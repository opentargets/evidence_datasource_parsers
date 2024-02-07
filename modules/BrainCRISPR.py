#!/usr/bin/env python
"""This module generates crispr evidence based on Brain Crispr resource."""
from __future__ import annotations

import gzip
import json
import logging
import re
from functools import reduce
from typing import TYPE_CHECKING, List
from urllib.parse import quote

from pyspark import SparkFiles
from pyspark.sql import functions as f
from pyspark.sql import types as t
from requests import post

if TYPE_CHECKING:
    from pyspark.sql import DataFrame, SparkSession


# create and configure main logger
logger = logging.getLogger(__name__)


class CRISPRBrain:
    """Parser for CRISPR Brain resource."""

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
        disease_mapping_url: str,
    ) -> None:
        """Initialize CRISPR Brain parser.

        Args:
            spark (SparkSession): Spark session object.
            disease_mapping_url (str): file containing the study to disease mapping
        """
        self.spark = spark

        # Adding disease mapping file to spark context:
        self.spark.sparkContext.addFile(disease_mapping_url)
        self.disease_mapping_file = disease_mapping_url.split("/")[-1]

    def __get_study_table(self: CRISPRBrain) -> DataFrame:
        """Extract screen data from Crispr brain API.

        Returns:
            DataFrame: with the parsed CRISPR Brain screens.

        Raises:
            ValueError: API versions don't match with the parser.
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
        """Read screen to disease and contrast mappings.

        Returns:
            DataFrame: with the parsed disease mappings.
        """
        disease_mappings_df = (
            self.spark.read.csv(
                SparkFiles.get(self.disease_mapping_file), sep="\t", header=True
            )
            # Splitting list of efos to list object:
            .withColumn(
                "diseaseFromSourceMappedId",
                f.split(f.col("diseaseFromSourceMappedId"), ", "),
            )
        )

        logger.info(
            "Number of studies with disease mapping: %s", disease_mappings_df.count()
        )
        return disease_mappings_df

    def __get_literature_mapping(self: CRISPRBrain) -> DataFrame:
        """Read screen publication to pmid mappings.

        Returns:
            DataFrame: with the parsed literature mappings.
        """
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
    def __parsing_experiment(description: str) -> dict[str, str | None]:
        """Parse free-text experiment description into structured data.

        Args:
            description (str): Free text screen description provided by the API.

        Returns:
            dict[str, str | None]: dictionary with 3 optional fields: "title", "experiment", and "analysis".
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
            .encode("ascii", "ignore")
            .decode("ascii")
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
                if description_lines[i] == "## Experiment":
                    experiment = description_lines[i + 1]
                if description_lines[i] == "## Analysis":
                    analysis = description_lines[i + 1]
        except TypeError:
            experiment = None
            analysis = None

        return {"title": title, "experiment": experiment, "analysis": analysis}

    @staticmethod
    def __QC_studies(df: DataFrame) -> DataFrame:
        """Report missing literature and disease annotation on the study table. Apply disease filter.

        Args:
            df (DataFrame): Study table

        Returns:
            DataFrame: Same schema, filtered for available disease mapping.
        """
        # Get the number of studies with no literature mapping:
        studies_wo_lit = df.filter(f.col("literature").isNull()).count()
        logger.info("Number of studies without literature mapping: %", studies_wo_lit)

        # Get the number of studies with no disease mapping:
        studies_wo_disease = df.filter(
            f.col("diseaseFromSourceMappedId").isNull()
        ).count()
        logger.info("Number of studies without disease mapping: %s", studies_wo_disease)

        # Get the number of studies with disease mapping:
        filtered_df = df.filter(f.col("diseaseFromSourceMappedId").isNotNull())
        logger.info("Number of studies with disease mapping: %s", filtered_df.count())

        # Return filtered dataframe:
        return filtered_df

    def __parse_gene_list(self: CRISPRBrain, study_id: str) -> DataFrame:
        """Parse gene list for a single screen and returns a dataframe with hit genes.

        The function filters out non-hit genes and renames some columns using the PySpark functions col and lit.

        Args:
            study_id (str): The study_id parameter is a string that represents the unique identifier for a
        single screen.

        Returns:
            DataFrame: hit genes from a single screen.
        """
        screen_url = f"https://storage.googleapis.com/crisprbrain.appspot.com/api-data/{quote(study_id)}.csv.gz"
        self.spark.sparkContext.addFile(screen_url)

        return (
            self.spark.read.csv(
                SparkFiles.get(f"{study_id}.csv.gz"), header=True, inferSchema=True
            )
            .filter(f.col("Hit Class") != "Non-Hit")
            .select(
                f.col("Gene").alias("targetFromSourceId"),
                f.col("P Value").alias("resourceScore"),
                f.when(f.col("Hit Class") == "Positive Hit", f.lit("upper tail"))
                .when(f.col("Hit Class") == "Negative Hit", f.lit("lower tail"))
                .alias("statisticalTestTail"),
                f.col("Phenotype").alias("log2FoldChangeValue"),
                f.lit(study_id).alias("studyId"),
            )
        )

    def __get_all_hit_genes(self, study_identifiers: List[str]) -> DataFrame:
        """Retrieve hit genes from all screens.

        Args:
            study_identifiers (List[str]): List of study identifiers.

        Returns:
            DataFrame: with hit genes from all screens.
        """
        # Parsing all gene lists for all screens and combine them into a single dataframe:
        return reduce(
            lambda df1, df2: df1.unionByName(df2),
            map(self.__parse_gene_list, study_identifiers),
        )

    def create_crispr_brain_evidence(self: CRISPRBrain) -> None:
        """Create evidence from CRISPR Brain resource."""
        # Reading
        crispr_studies = (
            # Reading studies:
            self.__get_study_table()
            # Joining with resolved literature:
            .join(self.__get_literature_mapping(), on="Reference Link", how="left")
            .select(
                # Adding constants:
                f.lit("affected_pathway").alias("datatypeId"),
                f.lit("crispr_screen").alias("datasourceId"),
                f.lit("crispr_brain").alias("projectId"),
                # Renaming fields:
                f.col("Screen Name").alias("studyId"),
                f.col("Libraries Screened").alias("crisprScreenLibrary"),
                f.col("studySummary.title").alias("studyOverview"),
                # Dropping experiment column for now:
                # f.col("studySummary.experiment").alias("experiment"),
                f.col("Cell Type").alias("cellType"),
                f.when(f.col("Genotype") != "WT", f.col("Genotype"))
                .otherwise(None)
                .alias("geneticBackground"),
                "Phenotype",
                "literature",
            )
            # Joining resolved literature:
            .join(self.__get_disease_mapping(), on="studyId", how="left")
            # Keep contrast if provided with the curation, otherwise keep the phenotype annotation:
            .withColumn("contrast", f.coalesce("contrast", "Phenotype"))
            # Do a QC and filter on the study table at this point:
            .transform(self.__QC_studies)
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
        """Get the evidence dataframe.

        Raises:
            Exception: AttributeError: Evidence not yet generated!

        Returns:
            DataFrame: Evidence dataframe.
        """
        try:
            return self.evidence
        except AttributeError as e:
            raise Exception("Evidence not yet generated!") from e

    def save_evidence(self: CRISPRBrain, filename: str) -> None:
        """Save evidence dataframe to file.

        Args:
            filename (str): File to save the evidence to.

        Raises:
            Exception: AttributeError: Evidence not yet generated!
        """
        format = filename.split(".")[-1]

        # If the filename doesn't contain a format, we default to parquet:
        if format == filename:
            format = "parquet"

        # Save data:
        try:
            self.evidence.write.mode("overwrite").format(format).save(filename)
        except AttributeError as e:
            raise Exception("Evidence not yet generated!") from e
