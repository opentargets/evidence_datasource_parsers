import logging
import os
from psutil import virtual_memory
import requests
import tempfile
import json
import gcsfs

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession, DataFrame
import pyspark.sql.functions as f
from pyspark.sql.types import StringType, StructField, StructType, ArrayType
from pyspark import SparkFiles
from psutil import virtual_memory


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

    def __init__(self, cell_passport_file: str, spark) -> None:
        self.cell_passport_file = cell_passport_file
        self.spark = spark
        self.cell_map = self.generate_map()

    def generate_map(self) -> None:
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
            * Tissue id is the UBERON identifier for the tissue, fetched from OLS
        """

        # loading cell line annotation data from Sanger:
        cell_df = (
            # The following option is required to correctly parse CSV records which contain newline characters.
            self.spark.read.option("multiline", True)
            .csv(self.cell_passport_file, header=True, sep=",", quote='"')
            .withColumn("biomarkerList", self.parse_msi_status(f.col("msi_status")))
            .select(
                f.col("model_name").alias("name"),
                f.col("model_id").alias("id"),
                f.col("tissue"),
                f.col("biomarkerList"),
            )
            .persist()
        )

        # Generating a unique set of tissues in a pandas dataframe:
        tissues = cell_df.select("tissue").distinct().toPandas()
        logging.info(f"Found {len(tissues)} tissues.")

        # Generating a unique set of cell lines in a pandas series:
        # TODO: UBERON IDs should not be mapped via REST queries. This should be done via joining with owl data.
        mapped_tissues = tissues.assign(
            tissueId=lambda df: df.tissue.apply(self.lookup_uberon)
        )
        logging.info(
            f"Found mapping for {len(mapped_tissues.loc[mapped_tissues.tissueId.notna()])} tissues."
        )

        # Converting to spark dataframe:
        mapped_tissues_spark = self.spark.createDataFrame(mapped_tissues)

        # Joining with cell lines:
        return (
            cell_df.join(mapped_tissues_spark, on="tissue", how="left")
            # Generating the diseaseCellLines object:
            .select(
                f.regexp_replace(f.col("name"), "-", "").alias("name"),
                "id",
                "biomarkerList",
                f.struct(["tissue", "name", "id", "tissueId"]).alias("diseaseCellLine"),
            )
        )

    def get_mapping(self):
        return self.cell_map

    @staticmethod
    def lookup_uberon(tissue_label: str) -> str:
        """Mapping tissue labels to tissue identifiers (UBERON codes) via the OLS API."""

        url = f"https://www.ebi.ac.uk/ols/api/search?q={tissue_label.lower()}&queryFields=label&ontology=uberon&exact=true"
        r = requests.get(url).json()

        if r["response"]["numFound"] == 0:
            return None
        else:
            return r["response"]["docs"][0]["short_form"]

    @staticmethod
    @f.udf(
        ArrayType(
            StructType(
                [
                    StructField("name", StringType(), nullable=False),
                    StructField("description", StringType(), nullable=False),
                ]
            )
        )
    )
    def parse_msi_status(status: str) -> list:
        """Based on the content of the MSI status, we generate the corresponding biomarker object."""

        if status == "MSI":
            return [{"name": "MSI", "description": "Microsatellite instable"}]
        if status == "MSS":
            return [{"name": "MSS", "description": "Microsatellite stable"}]
        else:
            return None


def read_path(path: str, spark_instance) -> DataFrame:
    """Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single TSV file; a single JSON file; a directory with JSON files; a directory with Parquet files."""
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
    if parquet_files:
        return spark_instance.read.parquet(path)


def import_trait_mappings() -> DataFrame:
    """Load the remote trait mappings file to a Spark dataframe."""

    remote_trait_mappings_url = "https://raw.githubusercontent.com/opentargets/curation/22.09.1/mappings/disease/manual_string.tsv"

    SparkSession.getActiveSession().sparkContext.addFile(remote_trait_mappings_url)

    return (
        SparkSession.getActiveSession()
        .read.csv(SparkFiles.get("manual_string.tsv"), header=True, sep="\t")
        .select(
            f.col("PROPERTY_VALUE").alias("diseaseFromSource"),
            f.element_at(f.split(f.col("SEMANTIC_TAG"), "/"), -1).alias(
                "diseaseFromSourceMappedId"
            ),
        )
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
