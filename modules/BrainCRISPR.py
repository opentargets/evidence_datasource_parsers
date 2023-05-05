#!/usr/bin/env python
"""This module generates crispr evidence based on Brain Crispr resource."""
from __future__ import annotations
from typing import TYPE_CHECKING

from urllib.parse import quote
from requests import post
import gzip
import json
import re
from functools import reduce
from typing import List

from pyspark.sql import types as t, functions as f
from pyspark import SparkFiles

if TYPE_CHECKING:
    from pyspark.sql import DataFrame, SparkSession


class CRISPRBrain:
    # CRISPR Brain API details. API key is required and also allows for checking for new release.
    API_SERVER_URL = "https://crisprbrain.org"
    API_CLIENT_ID = "92f263647c43525d3e4f181aa7e348f26e32129c0f827321d9261cc9765c56c0"
    API_VERSION = 1

    # Literature mappings:
    LITERATURE_MAPPING = [
        ("https://www.biorxiv.org/content/10.1101/2021.09.11.459904v1", ["ppr393667"]),
        ("https://pubmed.ncbi.nlm.nih.gov/30297964/", ["30297964"]),
        ("https://www.biorxiv.org/content/10.1101/2021.08.23.457400v1", ["PPR387167"]),
        ("https://doi.org/10.1101/2020.06.27.175679", ["34031600"]),
        ("https://pubmed.ncbi.nlm.nih.gov/31422865/", ["31422865"]),
        ("https://pubmed.ncbi.nlm.nih.gov/30449619/", ["30449619"]),
    ]

    def __init__(
        self: CRISPRBrain,
        spark: SparkSession,
        disease_mapping: str,
    ) -> None:
        """Initialize CRISPR Brain parser

        Args:
            spark (SparkSession):
            disease_mapping (str): file containing the study to disease mapping
        """
        self.spark = spark
        self.disease_mapping = disease_mapping

    def __get_study_table(self: CRISPRBrain) -> DataFrame:
        """Extract scree data from Crispr brain API.

        Returns:
            spark dataframe with the parsed CRISPR Brain screens.
        """

        # Establish connection to get study level metadata:
        response = post(
            f"{self.API_SERVER_URL}/api/screens", data={"client_id": self.API_CLIENT_ID}
        )

        # Retrieve study data from API:
        json_text = gzip.decompress(response.content).decode()
        screens = json.loads(json_text)

        # If the API version doesn't match, we raise error:
        if int(screens["__version"]) != self.API_VERSION:
            raise ValueError(
                f'API versions don\'t match: {self.API_VERSION} vs {screens["__version"]}. Parser needs to be updated.'
            )

        # Getting study level metadata:
        screen_data = [
            {**value, **value["metadata"]}
            for _, value in screens.items()
            if isinstance(value, dict)
        ]

        # The parsed study data is converted into a spark dataframe:
        return self.spark.createDataFrame(screen_data).withColumn(
            "studySummary", self.__parsing_experiment(f.col("Description"))
        )

    def __get_disease_mapping(self: CRISPRBrain) -> DataFrame:
        """Read screen to disease and contrast mappings."""
        return (
            self.spark.read.csv(self.disease_mapping, sep="\t", header=True)
            # Splitting list of efos to list object:
            .withColumn(
                "diseaseFromSourceMappedId",
                f.split(f.col("diseaseFromSourceMappedId"), ", "),
            )
        )

    def __get_literature_mapping(self: CRISPRBrain) -> DataFrame:
        """Read screen publication to pmid mappings."""
        return self.spark.createDataFrame(
            self.LITERATURE_MAPPING, ["Reference Link", "literature"]
        )

    @staticmethod
    @f.udf(
        t.StructType(
            [
                t.StructField("title", t.StringType()),
                t.StructField("experiment", t.StringType()),
                t.StructField("analysis", t.StringType()),
            ]
        )
    )
    def __parsing_experiment(description: str) -> dict:
        """Parse free-text experiment description into structured data.

        Args:
            description (str): Free text screen description provided by the API.

        Returns:
            dict: dictionary with 3 optional fields:
            - title
            - experiment
            - analysis
        """

        # Cleaning text:
        repl_patterns = [(r"\*+", ""), (r"\r", ""), (r"\t", ""), (r"\n+", "\n")]

        description = (
            reduce(
                lambda string, pattern: re.sub(pattern[0], pattern[1], string),
                repl_patterns,
                description,
            )
            # Dropping some nasty utf8 characters:
            .encode("ascii", "ignore").decode("ascii")
        )
        # Split
        description_lines = re.split(r"\n+", description.strip())
        # Title should be always available:
        title = (
            re.sub(r"#+\s+", "", description_lines[0])
            if "#" in description_lines[0]
            else description_lines[0]
        )

        experiment = None
        analysis = None

        try:
            for i in range(len(description_lines)):
                if "## Experiment" == description_lines[i]:
                    experiment = description_lines[i + 1]
                if "## Analysis" == description_lines[i]:
                    analysis = description_lines[i + 1]
        except TypeError:
            experiment = None
            analysis = None

        return {"title": title, "experiment": experiment, "analysis": analysis}

    def __parse_gene_list(self: CRISPRBrain, study_id: str) -> DataFrame:
        """
        This function parses a gene list for a single screen and returns a dataframe with hit genes.

        Args:
          study_id (str): The study_id parameter is a string that represents the unique identifier for a
        single screen.

        Returns:
          a DataFrame with hit genes from a single screen. The DataFrame contains columns
            - targetFromSource,
            - pValue,
            - statisticalTestTail,
            - resourceScore, and
            - studyId.

        The function filters out non-hit genes and renames some columns using the PySpark functions col and lit.
        """
        screen_url = f"https://storage.googleapis.com/crisprbrain.appspot.com/api-data/{quote(study_id)}.csv.gz"
        self.spark.sparkContext.addFile(screen_url)

        return (
            self.spark.read.csv(
                SparkFiles.get(f"{study_id}.csv.gz"), header=True, inferSchema=True
            )
            .filter(f.col("Hit Class") != "Non-Hit")
            .select(
                f.col("Gene").alias("targetFromSource"),
                f.col("P Value").alias("pValue"),
                f.when(f.col("Hit Class") == "Positive Hit", f.lit("upper tail"))
                .when(f.col("Hit Class") == "Negative Hit", f.lit("lower tail"))
                .alias("statisticalTestTail"),
                f.col("Phenotype").alias("resourceScore"),
                f.lit(study_id).alias("studyId"),
            )
        )

    def __get_all_hit_genes(self, study_identifiers: List[str]) -> DataFrame:
        """Retrieving hit genes from all screens."""

        # Parsing all gene lists for all screens and combine them into a single dataframe:
        return reduce(
            lambda df1, df2: df1.unionByName(df2),
            map(self.__parse_gene_list, study_identifiers),
        )

    def create_crispr_brain_evidence(self: CRISPRBrain) -> None:
        """"""
        # Reading
        crispr_studies = (
            # Reading studies:
            self.__get_study_table()
            # Joining with resolved literature:
            .join(self.__get_literature_mapping(), on="Reference Link", how="left")
            .select(
                # Adding constants:
                f.lit("sysbio").alias("datatypeId"),
                f.lit("crispr_screen").alias("datasourceId"),
                f.lit("crispr_brain").alias("projectId"),
                # Renaming fields:
                f.col("Screen Name").alias("studyId"),
                f.col("studySummary.title").alias("studyOverview"),
                f.col("studySummary.experiment").alias("experiment"),
                f.col("Cell Type").alias("cellType"),
                f.when(f.col("Genotype") != "WT", f.col("Genotype"))
                .otherwise(None)
                .alias("geneticBackground"),
                "Phenotype",
            )
            # Joining resolved literature:
            .join(self.__get_disease_mapping(), on="studyId", how="inner")
            # Keep contrast if provided with the curation, otherwise keep the phenotype annotation:
            .withColumn("contrast", f.coalesce("contrast", "Phenotype"))
            .drop("Phenotype")
            .persist()
        )

        # Extract all study identifiers into a list:
        study_identifiers = [
            x["studyId"] for x in crispr_studies.select("studyId").distinct().collect()
        ]

        # Based on the study identifiers, hit genes are retrieved:
        hit_genes = self.__get_all_hit_genes(study_identifiers)

        # Join genes and studies together + Explode diseases:
        self.evidence = (
            crispr_studies.join(hit_genes, on="studyId", how="left")
            .withColumn(
                "diseaseFromSourceMappedId",
                f.explode(f.col("diseaseFromSourceMappedId")),
            )
            .persist()
        )

    def get_dataframe(self: CRISPRBrain) -> DataFrame:
        try:
            return self.evidence
        except AttributeError as e:
            raise Exception("Evidence not yet generated!") from e

    def save_evidence(self: CRISPRBrain, filename: str) -> None:
        format = filename.split(".")[-1]

        if format == filename:
            format = "parquet"

        # Save data:
        self.evidence.write.mode("overwrite").format(format).save(filename)
