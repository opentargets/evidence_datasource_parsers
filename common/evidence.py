"""Shared functions for evidence generation."""
from __future__ import annotations

import json
import logging
import os
import sys
import tempfile
from typing import TYPE_CHECKING, Dict, Optional

import gcsfs
import pyspark.sql.functions as f
import yaml
from psutil import virtual_memory
from pyspark import SparkFiles
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession

if TYPE_CHECKING:
    from pyspark.sql import Column, DataFrame


def initialize_logger(
    name: str, log_file: Optional[str] = None, log_level: int = logging.INFO
) -> None:
    """Initialize the logger.

    Args:
        name (str): Name of the logger. This is typically the name of the module. Required to identify the logger.
        log_file (str): Path to the log file.
        log_level (int): log level eg. logging.INFO, logging.ERROR

    Returns:
        None
    """
    # Setting the format of the log messages:
    log_formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    # Setting up stream handler:
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(log_formatter)
    logger.addHandler(stream_handler)

    # If a log file is provided, add that handler too:
    if log_file is not None:
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(log_formatter)
        logger.addHandler(file_handler)


def detect_spark_memory_limit():
    """Spark does not automatically use all available memory on a machine. When working on large datasets, this may
    cause Java heap space errors, even though there is plenty of RAM available. To fix this, we detect the total amount
    of physical memory and allow Spark to use (almost) all of it."""
    mem_gib = virtual_memory().total >> 30
    return int(mem_gib * 0.9)


def write_evidence_strings(evidence, output_file):
    """Exports the table to a compressed JSON file containing the evidence strings."""
    with tempfile.TemporaryDirectory() as tmp_dir_name:
        (
            evidence.coalesce(1)
            .write.format("json")
            .mode("overwrite")
            .option("compression", "org.apache.hadoop.io.compress.GzipCodec")
            .save(tmp_dir_name)
        )
        json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith(".json.gz")]
        assert (
            len(json_chunks) == 1
        ), f"Expected one JSON file, but found {len(json_chunks)}."
        os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)


def initialize_sparksession() -> SparkSession:
    """Initialize spark session."""

    spark_mem_limit = detect_spark_memory_limit()
    spark_conf = (
        SparkConf()
        .set("spark.driver.memory", f"{spark_mem_limit}g")
        .set("spark.executor.memory", f"{spark_mem_limit}g")
        .set("spark.driver.maxResultSize", "0")
        .set("spark.debug.maxToStringFields", "2000")
        .set("spark.sql.execution.arrow.maxRecordsPerBatch", "500000")
        .set("spark.ui.showConsoleProgress", "false")
    )
    return (
        SparkSession.builder.config(conf=spark_conf)
        .master("local[*]")
        .config("spark.driver.bindAddress", "127.0.0.1")
        .getOrCreate()
    )


class GenerateDiseaseCellLines:
    """
    Generate "diseaseCellLines" object from a cell passport file.

    !!!
    There's one important bit here: I have noticed that we frequenty get cell line names
    with missing dashes. Therefore the cell line names are cleaned up by removing dashes.
    It has to be done when joining with other datasets.
    !!!

    Args:
        cell_passport_file: Path to the cell passport file.
    """

    def __init__(
        self,
        spark: SparkSession,
        cell_passport_file: str,
        cell_line_to_uberon_mapping: str,
    ) -> None:
        self.cell_passport_file = cell_passport_file
        self.spark = spark
        self.tissue_to_uberon_map = self.read_uberon_mapping(
            cell_line_to_uberon_mapping
        )

    def generate_disease_cell_lines(self) -> DataFrame:
        """Reading and procesing cell line data from the cell passport file.

        The schema of the returned dataframe is:

        root
        |-- name: string (nullable = true)
        |-- id: string (nullable = true)
        |-- biomarkerList: array (nullable = true)
        |    |-- element: struct (containsNull = true)
        |    |    |-- name: string (nullable = true)
        |    |    |-- description: string (nullable = true)
        |-- diseaseCellLine: struct (nullable = false)
        |    |-- tissue: string (nullable = true)
        |    |-- name: string (nullable = true)
        |    |-- id: string (nullable = true)
        |    |-- tissueId: string (nullable = true)

        Note:
            * Microsatellite stability is the only inferred biomarker.
            * The cell line name has the dashes removed.
            * Id is the cell line identifier from Sanger
            * Tissue id is the UBERON identifier for the tissue, based on manual curation.
        """
        # loading cell line annotation data from Sanger:
        cell_df = (
            # The following option is required to correctly parse CSV records which contain newline characters.
            self.spark.read.option("multiline", True)
            .csv(self.cell_passport_file, header=True, sep=",", quote='"')
            .select(
                f.col("model_name").alias("name"),
                f.col("model_id").alias("id"),
                f.lower(f.col("tissue")).alias("tissueFromSource"),
                f.array(self.parse_msi_status(f.col("msi_status"))).alias(
                    "biomarkerList"
                ),
            )
            # Joning with the UBERON mapping:
            .join(self.tissue_to_uberon_map, on="tissueFromSource", how="left")
            .persist()
        )

        # Joining with cell lines:
        return (
            cell_df
            # Generating the diseaseCellLines object:
            .select(
                f.regexp_replace(f.col("name"), "-", "").alias("name"),
                "id",
                "biomarkerList",
                f.struct(
                    f.col("tissueName").alias("tissue"), "name", "id", "tissueId"
                ).alias("diseaseCellLine"),
            )
        )

    def read_uberon_mapping(self, cell_line_to_uberon_mapping: str) -> DataFrame:
        """Read the cell line to UBERON mapping file into a Spark dataframe (http or local).

        Returned schema:
            root
            |-- tissue: string (nullable = true)
            |-- tissueId: string (nullable = true)
            |-- tissueName: string (nullable = true)

        Args:
            cell_line_to_uberon_mapping (str): Path to the cell line to UBERON mapping file.

        Returns:
            DataFrame: Spark dataframe containing the cell line to UBERON mapping.
        """

        # If the mapping file is a URL, download it and read it into a Spark dataframe.
        if "http" in cell_line_to_uberon_mapping:
            self.spark.sparkContext.addFile(cell_line_to_uberon_mapping)
            cell_line_to_uberon_mapping = SparkFiles.get(
                cell_line_to_uberon_mapping.split("/")[-1]
            )

        # Reading the mapping file into a Spark dataframe.
        return self.spark.read.csv(cell_line_to_uberon_mapping, header=True, sep=",")

    @staticmethod
    def parse_msi_status(status: Column) -> Column:
        """Based on the content of the MSI status, we generate the corresponding biomarker object."""

        return (
            f.when(
                status == "MSI",
                f.struct(
                    f.lit("MSI").alias("name"),
                    f.lit("Microsatellite instable").alias("description"),
                ),
            )
            .when(
                status == "MSS",
                f.struct(
                    f.lit("MSS").alias("name"),
                    f.lit("Microsatellite stable").alias("description"),
                ),
            )
            .otherwise(f.lit(None))
        )


def read_path(path: str, spark_instance) -> DataFrame:
    """Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single TSV file; a single JSON file; a directory with JSON files; a directory with Parquet files.
    """
    if path is None:
        return None

    # The provided path must exist and must be either a file or a directory.
    assert os.path.exists(path), f"The provided path {path} does not exist."
    assert os.path.isdir(path) or os.path.isfile(
        path
    ), f"The provided path {path} is neither a file or a directory."

    # Case 1: We are provided with a single file.
    if os.path.isfile(path):
        if path.endswith(".csv"):
            return spark_instance.read.csv(path, header=True, inferSchema=True)
        if path.endswith(".tsv"):
            return spark_instance.read.csv(path, sep="\t", header=True)
        elif path.endswith((".json", ".json.gz", ".jsonl", ".jsonl.gz")):
            return spark_instance.read.json(path)
        else:
            raise AssertionError(
                f"The format of the provided file {path} is not supported."
            )

    # Case 2: We are provided with a directory. Let's peek inside to see what it contains.
    all_files = [
        os.path.join(dp, filename)
        for dp, dn, filenames in os.walk(path)
        for filename in filenames
    ]

    # It must be either exclusively JSON, or exclusively Parquet.
    json_files = [
        fn
        for fn in all_files
        if fn.endswith((".json", ".json.gz", ".jsonl", ".jsonl.gz"))
    ]
    parquet_files = [fn for fn in all_files if fn.endswith(".parquet")]
    assert not (
        json_files and parquet_files
    ), f"The provided directory {path} contains a mix of JSON and Parquet."
    assert (
        json_files or parquet_files
    ), f"The provided directory {path} contains neither JSON nor Parquet."

    # A directory with JSON files.
    if json_files:
        return spark_instance.read.option("recursiveFileLookup", "true").json(path)

    # A directory with Parquet files.
    else:
        return spark_instance.read.parquet(path)


def read_project_config() -> Dict[str, str]:
    """Load the configuration file into a dictionary."""
    base_dir = os.getcwd()
    if "common" in base_dir:
        config_path = os.path.join(os.path.dirname(base_dir), "configuration.yaml")
    else:
        config_path = os.path.join(base_dir, "configuration.yaml")
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file ({config_path}) does not exist.")
    return yaml.safe_load(open(config_path))


def import_trait_mappings(spark: SparkSession) -> DataFrame:
    """Load the remote trait mappings file to a Spark dataframe. The file is downloaded from the curated data repository.

    Args:
        spark (SparkSession): Spark session.

    Returns:
        DataFrame: Spark dataframe containing `diseaseFromSource` and `diseaseFromSourceMappedId` columns.
    """
    remote_trait_mappings_url = f"{read_project_config()['global']['curation_repo']}/mappings/disease/manual_string.tsv"
    spark.sparkContext.addFile(remote_trait_mappings_url)
    return spark.read.csv(
        SparkFiles.get("manual_string.tsv"), header=True, sep="\t"
    ).select(
        f.col("PROPERTY_VALUE").alias("diseaseFromSource"),
        f.element_at(f.split(f.col("SEMANTIC_TAG"), "/"), -1).alias(
            "diseaseFromSourceMappedId"
        ),
    )


def read_ppp_config(config_path: str) -> dict:
    """Read json file stored on GCP location into a dictionary.

    Args:
        config_path (str): Path to configuration file.

    Returns:
        dict: parsed configuration
    """
    # Configuration is provided in a gs:// location:
    if config_path.startswith("gs://"):
        fs = gcsfs.GCSFileSystem()
        try:
            with fs.open(config_path, "r") as parameter_file:
                parameters = json.load(parameter_file)
        except Exception as e:
            raise e(f"Could not read parameter file. {config_path}")

    # If not a GCP location, assuming a local file:
    else:
        try:
            with open(config_path, "r") as parameter_file:
                parameters = json.load(parameter_file)
        except Exception as e:
            raise e(f"Could not read parameter file. {config_path}")

    return parameters
