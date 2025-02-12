#!/usr/bin/env python3
"""Parser to generate evidence for ot_crispr datasets.

This dataset consist of a series of OTAR projects studying various diseases using genome-wide crisp/cas9 knock-outs.
 - The results are expected to arrive in MAGeCK format.
 - The study level metadata is expected to come via filling out a Google spreadseet.
 - These spreadseet downloaded as a tsv and version-ed in the PPP-evidencie-configuration repository.
 - The generated evidence is exptected to be validated against the OT evidence schema.
"""

from __future__ import annotations

import argparse
import logging
from functools import reduce
from typing import Optional

from pyspark.sql import Column, DataFrame, Row, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t

from src.common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)


class OTAR_CRISPR_study_parser:
    """Process raw CRISPR study table.

    - The versioned, raw study table provided as tsv.
    - The following operations are performed:
        - Read the study table.
        - Drop rows with non-OTAR projects.
        - Split columns with multiple values for disease, filterColumn, and dataFile.
        - Collect replicates for each study.
    """

    def __init__(self, spark: SparkSession, study_table_path: str) -> None:
        """Initialise the study parser.

        Args:
            spark (SparkSession): _description_
            study_table_path (str): path to the study tsv.
        """
        self.study_table_path = study_table_path
        self.spark = spark
        self.logger = logging.getLogger(__name__)
        self.study_table: DataFrame

    @staticmethod
    def split_column_value(col: Column, separator: str = r"\|") -> Column:
        """Remove whitespace and split column value by the provided separator.

        Args:
            col (Column): A column to be split.
            separator (str): A separator to split the column value.

        Returns:
            Column: A column with split values.
        """
        return f.split(f.regexp_replace(col, r"^\s+", ""), separator)

    def log_study_table(self: OTAR_CRISPR_study_parser) -> None:
        """Report key metrics on study table.

        Args:
            self (OTAR_CRISPR_study_parser).
        """
        # Get study count:
        self.logger.info(f"Number of studies processed: {self.study_table.count()}")
        # Get project count:
        self.logger.info(
            f'Number of projects: {self.study_table.select("projectId").distinct().count()}'
        )

        # Get studies with multiple replicates:
        studies_with_multiple_replicates = self.study_table.filter(
            f.size("replicates") > 1
        ).count()
        self.logger.info(
            f"Studies with multiple replicates: {studies_with_multiple_replicates}"
        )

    def process_studies(self: OTAR_CRISPR_study_parser) -> OTAR_CRISPR_study_parser:
        """Parsing the study table for ot_crispr datasets.

        Args:
            study_table (DataFrame): A DataFrame with study level metadata.

        Returns:
            DataFrame: A DataFrame with parsed study level metadata.
        """
        study_table = (
            self.spark.read.csv(
                self.study_table_path, header=True, inferSchema=True, sep="\t"
            )
            # Dropping studies with no OTAR project and the field description row:
            .filter(f.col("projectId").startswith("OTAR"))
            # Selecting relevant columns:
            .select(
                "studyId",
                "projectId",
                "projectDescription",
                "studyOverview",
                "releaseVersion",
                f.col("releaseDate").cast(t.DateType()).alias("releaseDate"),
                # Splitting diseases:
                self.split_column_value(f.col("diseases")).alias(
                    "diseaseFromSourceMappedId"
                ),
                "isCellTypeDerived",
                "crisprScreenLibrary",
                "crisprStudyMode",
                "geneticBackground",
                "cellType",
                "cellLineBackground",
                "contrast",
                "dataFileType",
                # Splitting and exploding study if both tail of the distribution are used:
                f.explode_outer(
                    self.split_column_value(f.col("filterColumn"), ",")
                ).alias("filterColumn"),
                # Casting threshold to float:
                f.col("threshold").cast(t.FloatType()).alias("threshold"),
                "dataFile",
                "ControlDataset",
                # Adding replicate identifier when missing:
                f.when(f.col("replicateNumber").isNull(), f.lit(1))
                .otherwise(f.col("replicateNumber"))
                .alias("replicateId"),
            )
            # Grouping by study level:
            .groupBy(
                "studyId",
                "projectId",
                "projectDescription",
                "studyOverview",
                "releaseVersion",
                "releaseDate",
                "diseaseFromSourceMappedId",
                "isCellTypeDerived",
                "crisprScreenLibrary",
                "crisprStudyMode",
                "geneticBackground",
                "cellType",
                "cellLineBackground",
                "contrast",
                "filterColumn",
                "threshold",
            )
            # Collecting replicates for each study:
            .agg(
                f.collect_list(
                    f.struct("dataFile", "ControlDataset", "replicateId")
                ).alias("replicates")
            )
        )
        assert isinstance(study_table, DataFrame), "Failed to parse the study table."
        self.study_table = study_table

        # Let's check the study table:
        self.log_study_table()
        return self

    def get_study_data(self: OTAR_CRISPR_study_parser) -> DataFrame:
        """Return the study table."""
        return self.study_table


class OTAR_CRISPR_evidence_generator:
    """Generate evidence from OTAR CRISPR datasets.

    Based on the provided study-level metadata and the provided raw data files, the evidence is generated.

    The following operations are performed:
        - Process replicates for each study.
        - Combine the results into a single DataFrame.
        - Add study level metadata.
        - Write the evidence to a file.
    """

    DATASOURCE_ID = "ot_crispr"
    DATATYPE_ID = "ot_partner"

    def __init__(
        self, spark: SparkSession, study_table: DataFrame, data_path: str
    ) -> None:
        """Initialise the evidence generator."""
        self.spark = spark
        self.study_table = study_table
        self.logger = logging.getLogger(__name__)
        self.data_path = data_path

    @staticmethod
    def _get_test_side(filter_column: str) -> Column:
        """Get the test side based on the filter column.

        Args:
            filter_column (str): A filter column name.

        Returns:
            Column: A column with the test side.
        """
        return f.when(f.lit(filter_column).contains("pos"), f.lit("upper tail")).when(
            f.lit(filter_column).contains("neg"), f.lit("lower tail")
        )

    def _read_and_filter_mageck_file(
        self: OTAR_CRISPR_evidence_generator,
        mageck_file: str,
        filter_column: str,
        threshold: float,
    ) -> DataFrame:
        """Read and filter MAGeCK files based on the provided threshold applied on the specified column.

        Args:
            mageck_file (str): A list of files to be read.
            filter_column (str): A filter column name.
            threshold (float): A threshold to filter the data.

        Returns:
            DataFrame: A DataFrame with filtered data.

        Raises:
            ValueError: If the label separator is not recognized.
        """
        side = filter_column.split("|")[0]

        # Some MAGEcK files have different column label separators eg. "pos|lfc" vs "pos.lfc" We have to sort this out:
        label_separator = "|"
        raw_data = self.spark.read.csv(mageck_file, header=True, sep="\t")

        # Checking label separator in the third column, which expected to be: neg|p-value or neg.p-value:
        if "|" in raw_data.columns[3]:
            label_separator = "|"
        elif "." in raw_data.columns[3]:
            label_separator = "."
        else:
            raise ValueError(f"Unrecognized label separator in {raw_data.columns[2]}")

        # Updating column names according to the identified label separator:
        raw_data = reduce(
            # Rename all columns:
            lambda df, col: df.withColumnRenamed(
                col, col.replace(label_separator, "|")
            ),
            raw_data.columns,
            raw_data,
        )

        return raw_data.select(
            f.split(f.col("id"), "_")[0].alias("targetFromSourceId"),
            f.col(f"{side}|lfc").cast(t.FloatType()).alias("log2FoldChangeValue"),
            f.col(filter_column).cast(t.FloatType()).alias("resourceScore"),
            self._get_test_side(filter_column).alias("statisticalTestTail"),
        ).filter(f.col("resourceScore") < threshold)

    def _process_replicate(
        self: OTAR_CRISPR_evidence_generator,
        data_file: str,
        control_dataset: Optional[str],
        filter_column: str,
        threshold: float,
    ) -> DataFrame:
        """Process a single replicate: finding hits, exclude controls if provided.

        Args:
            data_file (str): A single file in mageck format.
            control_dataset (Optional[str]): A control dataset in mageck format.
            filter_column (str): A filter column name.
            threshold (float): A threshold to filter the data.

        Returns:
            DataFrame: A DataFrame with processed data.
        """
        # Extract hist from the data files:
        hits = self._read_and_filter_mageck_file(data_file, filter_column, threshold)
        # If control dataset is provided, filter out the hits:
        if control_dataset:
            hits = hits.join(
                (
                    self._read_and_filter_mageck_file(
                        control_dataset, filter_column, threshold
                    )
                    .select("targetFromSourceId")
                    .distinct()
                ),
                how="left_anti",
                on="targetFromSourceId",
            )
        return hits

    def _process_study_table_row(
        self: OTAR_CRISPR_evidence_generator, row: Row
    ) -> DataFrame:
        """Process a single row from the study table.

        Args:
            row (Row): A row from the study table.

        Returns:
            DataFrame: A DataFrame with processed data.
        """
        # Process all replicates and collect the results in a list of dataframes:
        replicate_data = [
            self._process_replicate(
                data_file=f"{self.data_path}/{row.projectId}/{replicate.dataFile}",
                control_dataset=f"{self.data_path}/{row.projectId}/{replicate.ControlDataset}"
                if replicate.ControlDataset
                else None,
                filter_column=row.filterColumn,
                threshold=row.threshold,
            )
            for replicate in row.replicates
        ]

        return (
            # Combine the results into a single DataFrame:
            reduce(lambda df1, df2: df1.unionByName(df2), replicate_data)
            # Aggregating replicate level data:
            .groupBy("targetFromSourceId")
            .agg(
                f.collect_list(
                    f.struct(
                        "log2FoldChangeValue", "resourceScore", "statisticalTestTail"
                    )
                ).alias("replicates")
            )
            # Add replicate count:
            .withColumn("replicateCount", f.lit(len(replicate_data)))
            # Drop genes, which were not found in all replicates:
            .filter(f.size("replicates") == f.col("replicateCount"))
            # Select the best replicate:
            .select(
                "targetFromSourceId",
                f.col("replicates")[0].log2FoldChangeValue.alias("log2FoldChangeValue"),
                f.col("replicates")[0].resourceScore.alias("resourceScore"),
                f.col("replicates")[0].statisticalTestTail.alias("statisticalTestTail"),
                f.lit(row.studyId).alias("studyId"),
            )
        )

    def save_evidence_data(
        self: OTAR_CRISPR_evidence_generator, output_file: str
    ) -> None:
        # Process the study table:
        all_hits = reduce(
            lambda df1, df2: df1.unionByName(df2),
            [
                self._process_study_table_row(study)
                for study in self.study_table.collect()
            ],
        )

        # Adding study level metadata:
        evidence = (
            all_hits
            # Joining study level metadata:
            .join(self.study_table, on="studyId", how="inner")
            # Selecting relevant columns:
            .select(
                # Project level metadata:
                "projectId",
                "projectDescription",
                "releaseDate",
                "releaseVersion",
                # Study level metadata:
                "studyId",
                "studyOverview",
                "contrast",
                "crisprScreenLibrary",
                "cellType",
                "cellLineBackground",
                f.explode("diseaseFromSourceMappedId").alias(
                    "diseaseFromSourceMappedId"
                ),
                # Evidence level data:
                "targetFromSourceId",
                "log2FoldChangeValue",
                "resourceScore",
                "statisticalTestTail",
                # Static fields:
                f.lit("ot_crispr").alias("datasourceId"),
                f.lit("ot_partner").alias("datatypeId"),
            )
        )
        self.logger.info(f"Number of evidence: {evidence.count()}")
        write_evidence_strings(evidence, output_file)
        self.logger.info("Done.")


def main(study_table_path, output_file, data_folder) -> None:
    # Get logger:
    logger = logging.getLogger(__name__)

    # Initialize spark session:
    spark = initialize_sparksession()

    # Report input data:
    logger.info(f"Study information is read from: {study_table_path}")
    logger.info(f"Evidence saved to: {output_file}")
    logger.info(f"Data files read from: {data_folder}")

    # Parsing study table:
    study_table = (
        OTAR_CRISPR_study_parser(spark, study_table_path)
        .process_studies()
        .get_study_data()
    )

    # Parsing study table:
    OTAR_CRISPR_evidence_generator(spark, study_table, data_folder).save_evidence_data(
        output_file
    )


def parse_arguments() -> argparse.Namespace:
    """Parsing command line arguments.

    Returns:
        argparse.Namespace: Parser for the ot_crispr evidence generation
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--study_table",
        type=str,
        required=True,
        help="A JSON file with study level metadata.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output json.gz file with the evidece.",
    )
    parser.add_argument(
        "-l",
        "--log_file",
        required=False,
        type=str,
        help="Logs are saved into this file.",
    )
    parser.add_argument(
        "-d",
        "--data_folder",
        required=True,
        type=str,
        help="Folder with the data files.",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Parse arguments:
    args = parse_arguments()

    study_table = args.study_table
    log_file = args.log_file
    data_folder = args.data_folder
    output_file = args.output

    # Configure logger:
    initialize_logger(__name__, log_file)

    # Process data:
    main(study_table, output_file, data_folder)
