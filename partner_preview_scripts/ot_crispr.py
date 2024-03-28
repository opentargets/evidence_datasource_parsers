#!/usr/bin/env python3
"""Parser to generate evidence for ot_crispr datasets.

This dataset consist of a series of OTAR projects studying various diseases using genome-wide crisp/cas9 knock-outs.
 - The results are expected to arrive in MAGeCK format.
 - The study level metadata is expected to come via filling out a Google spreadseet.
 - These spreadseet is then converted into a json and version-ed in the PPP-evidencie-configuration repository.
 - The generated evidence is exptected to be validated against the OT evidence schema.

"""

from __future__ import annotations

import logging
import sys
from datetime import datetime
from functools import reduce
from typing import List, Optional

import pandas as pd
from pyspark.sql import Column, DataFrame, Row, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)


class OTAR_CRISPR_study_parser_old:
    def __init__(self, study_url: str) -> None:
        self.study_file = study_url
        self.logger = logging.getLogger(__name__)

        # Store the study dataframe after dropping problematic studies:
        study_df = (
            pd.read_json(study_url)
            # drop rows with no study id or data file:
            .loc[lambda df: df.studyId.notna() & df.dataFiles.notna()]
        )

        # Test and warn if multiple studies have the same study id/replicateNumber pairs:
        duplicated_study_ids = (
            study_df.assign(
                study_replicate=lambda df: df.studyId
                + df.replicateNumber.fillna("").astype(str)
            )
            .loc[lambda df: df.study_replicate.duplicated()]
            .studyId.tolist()
        )

        assert (
            len(duplicated_study_ids) == 0
        ), f'Multiple studies have the same study id without distinguishing replicate identifier: {", ".join(duplicated_study_ids)}'

        # Get the number of replicates for each study:
        replicates = (
            study_df.groupby("studyId")
            .agg(replicateCount=pd.NamedAgg(column="studyId", aggfunc="count"))
            .reset_index()
        )

        # Splitting studies if
        study_df = (
            study_df
            # Joining with replicate number:
            .merge(replicates, on="studyId", how="left")
            # Exploding filter columns:
            .assign(filterColumn=lambda df: df.filterColumn.str.split(","))
            .explode("filterColumn")
        )

        self.logger.info(
            f"Number of studies processed: {len(study_df.studyId.unique())}"
        )

        projects = study_df.projectId.unique()
        self.logger.info(f'Number of projects: {len(projects)} ({", ".join(projects)})')

        self.study_df = study_df

    def generate_evidence(self, data_folder: str) -> None:
        # Looping through the studies and generating evidence:
        # Reading all data files and filter for significant hits:
        study_columns = [
            "releaseDate",
            "releaseVersion",
            "studyId",
            "dataFiles",
            "dataFileType",
            "filterColumn",
            "threshold",
            "projectId",
            "ControlDataset",
            "projectDescription",
            "replicateNumber",
        ]

        # hits is a pd.Series with pd.DataFrames as values.
        hits = (
            self.study_df.explode("dataFiles")
            .assign(
                dataFile=lambda df: df.apply(
                    lambda x: f'{data_folder}/{x["projectId"]}/{x["dataFiles"]}', axis=1
                )
            )
            .assign(
                ControlDataset=lambda df: df.apply(
                    lambda x: f'{data_folder}/{x["projectId"]}/{x["ControlDataset"]}'
                    if pd.notna(x["ControlDataset"])
                    else None,
                    axis=1,
                )
            )
            # TODO: parsing the data files should be file type dependent!
            # The following apply returns pd.DataFrames:
            .apply(self.parse_MAGeCK_file, axis=1)
        )

        # Concatenate all hits into one single dataframe:
        hits_df = (
            pd.concat(hits.to_list())
            .reset_index(drop=True)
            # Cleaning gene id column:
            .assign(
                targetFromSourceId=lambda df: df.targetFromSourceId.apply(
                    self.cleaning_gene_id
                )
            )
        )

        # Merging:
        evidence_fields = [
            "targetFromSourceId",
            "diseaseFromSourceMappedId",
            "projectDescription",
            "projectId",
            "studyId",
            "studyOverview",
            "contrast",
            "crisprScreenLibrary",
            "cellType",
            "cellLineBackground",
            "geneticBackground",
            "statisticalTestTail",
            "resourceScore",
            "log2FoldChangeValue",
            "releaseDate",
            "releaseVersion",
        ]
        self.merged_dataset = (
            self.study_df.assign(
                statisticalTestTail=lambda df: df.filterColumn.map(FILTER_COLUMN_MAP)
            )
            .merge(hits_df, on=["studyId", "filterColumn"], how="inner")
            .explode("diseaseFromSourceMappedId")
            # .filter(items=evidence_fields)
            .assign(datasourceId="ot_crispr", datatypeId="ot_partner")
        )

        # Save the merged dataset to file:
        self.merged_dataset.to_csv("merged_dataset.csv")

    @staticmethod
    def cleaning_gene_id(gene_id: str) -> str:
        """Expandable set of string processing steps to clean gene identifiers.

        Examples:
            >>> cleaning_gene_id("ENSG00000187123_LYPD6")
            >>> "ENSG00000187123"
        """

        # ENSG00000187123_LYPD6 -> ENSG00000187123
        gene_id = gene_id.split("_")[0]

        return gene_id

    @staticmethod
    def parse_MAGeCK_file(row: pd.Series) -> pd.DataFrame:
        """This function returns a pandas dataframe with the datafile and with properly named columns"""

        datafile = row["dataFile"]
        filterColumn = row["filterColumn"]
        threshold = float(row["threshold"])
        studyId = row["studyId"]
        controlDataFile = row["ControlDataset"]

        # Which end of the distribution are we looking? - "neg" or "pos"?
        side = filterColumn.split("|")[0]

        # Read data, filter and rename columns:
        mageck_df = (
            pd.read_csv(datafile, sep="\t")
            .rename(
                columns={
                    filterColumn: "resourceScore",
                    "id": "targetFromSourceId",
                    # Extracting log fold change for the relevant direction:
                    f"{side}|lfc": "log2FoldChangeValue",
                }
            )
            .loc[lambda df: df.resourceScore <= threshold][
                ["targetFromSourceId", "resourceScore", "log2FoldChangeValue"]
            ]
            .assign(studyId=studyId, filterColumn=filterColumn)
        )

        # Applying control if present:
        if pd.isna(controlDataFile):
            logging.info(f"Number of genes reach threshold: {len(mageck_df)}")
            return mageck_df

        # Read control data, filter and rename columns:
        logging.info(f"Reading control data file: {controlDataFile}")
        controlHits = (
            pd.read_csv(controlDataFile, sep="\t")
            .rename(columns={filterColumn: "resourceScore", "id": "targetFromSourceId"})
            .loc[lambda df: df.resourceScore <= threshold]["targetFromSourceId"]
            .tolist()
        )

        # Excluding control genes:
        mageck_df = mageck_df.loc[lambda df: df.targetFromSourceId.isin(controlHits)]

        logging.info(f"Number of genes reach threshold: {len(mageck_df)}")
        return mageck_df

    def write_evidence(self, output_file: str) -> None:
        """Write the merged evidence to file"""

        json_list = [
            json.dumps(row.dropna().to_dict())
            for _, row in self.merged_dataset.iterrows()
        ]
        with gzip.open(output_file, "wt") as f:
            f.write("\n".join(json_list) + "\n")


class OTAR_CRISPR_study_parser:
    def __init__(self, spark: SparkSession, study_table_path: str) -> None:
        self.study_table_path = study_table_path
        self.spark = spark
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def split_column_value(col: Column, separator: str = r"\|") -> Column:
        """Remove whitespace and split column value by separator.

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
        studies_with_multiple_replicates = study_table.filter(
            f.size("replicates") > 1
        ).count()
        self.logger.info(
            f"Studies with multiple replicates: {studies_with_multiple_replicates}"
        )

    def process_studies(self: OTAR_CRISPR_study_parser) -> None:
        """Parsing the study table for ot_crispr datasets.

        Args:
            study_table (DataFrame): A DataFrame with study level metadata.

        Returns:
            DataFrame: A DataFrame with parsed study level metadata.
        """
        # Helper function to split column values:

        self.study_table = (
            self.spark.read.csv(
                self.study_table_path, header=True, inferSchema=True, sep="\t"
            )
            # Dropping rows not starting with OTAR project:
            .filter(f.col("projectId").startswith("OTAR"))
            # Filling replicateId:
            .select(
                "studyId",
                "projectId",
                "projectDescription",
                "studyOverview",
                "releaseVersion",
                "releaseDate",
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
                f.explode_outer(
                    self.split_column_value(f.col("filterColumn"), ",")
                ).alias("filterColumn"),
                f.col("threshold").cast(t.FloatType()).alias("threshold"),
                self.split_column_value(f.col("dataFile")).alias("dataFiles"),
                "ControlDataset",
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
                    f.struct("dataFiles", "ControlDataset", "replicateId")
                ).alias("replicates")
            )
        )
        # Let's check the study table:
        self.log_study_table()
        return self

    def get_study_data(self: OTAR_CRISPR_study_parser) -> DataFrame:
        """Return the study table."""
        return self.study_table


class OTAR_CRISPR_evience_generator:
    DATASOURCE_ID = "ot_crispr"
    DATATYPE_ID = "ot_partner"

    def __init__(
        self, spark: SparkSession, study_table: DataFrame, data_path: str
    ) -> None:
        self.spark = spark
        self.study_table = study_table
        self.logger = logging.getLogger(__name__)
        self.data_path = data_path
        self.data_columns = [
            "studyId",
            "projectId",
            "filterColumn",
            "threshold",
            "replicates",
        ]

    @staticmethod
    def _get_test_side(filter_column: str) -> Column:
        """Get the test side based on the filter column.

        Args:
            filter_column (str): A filter column name.

        Returns:
            Column: A column with the test side.
        """
        return f.when(f.col("filterColumn").contains("pos"), f.lit("upper tail")).when(
            f.col("filterColumn").contains("neg"), f.lit("lower tail")
        )

    def _read_and_filter_mageck_files(
        self: OTAR_CRISPR_evience_generator,
        files: List[str],
        filter_column: str,
        threshold: float,
    ) -> DataFrame:
        side = filter_column.split("|")[0]
        return (
            self.spark.read.csv(files, header=True, sep="\t")
            .select(
                f.col("id").alias("targetFromSourceId"),
                f.col(f"{side}|lfc").alias("log2FoldChangeValue"),
                f.col(filter_column).alias("resourceScore"),
                self._get_test_side(filter_column).alias("statisticalTestTail"),
            )
            .filter(f.col("resourceScore") < threshold)
            .persist()
        )

    def _process_replicate(
        self: OTAR_CRISPR_evience_generator,
        data_files: List[str],
        control_dataset: Optional[str],
        filter_column: str,
        threshold: float,
    ):
        hits = self._read_and_filter_mageck_files(data_files, filter_column, threshold)
        if control_dataset:
            hits = hits.join(
                (
                    self._read_and_filter_mageck_files(
                        [control_dataset], filter_column, threshold
                    )
                    .select("targetFromSourceId")
                    .distinct()
                ),
                how="left_anti",
                on="targetFromSourceId",
            )
        return hits

    def _process_process_study_table_row(
        self: OTAR_CRISPR_evience_generator, row: Row
    ) -> DataFrame:
        # Process all replicates and collect the results in a list of dataframes:
        replicate_data = [
            self._process_replicate(
                data_files=[
                    f"{self.data_path}/{row.projectId}/{data_file}"
                    for data_file in replicate.dataFiles
                ],
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

    def generate_evidence(self: OTAR_CRISPR_evience_generator) -> Optional[DataFrame]:
        # Process the study table:
        all_hits = reduce(
            lambda df1, df2: df1.unionByName(df2),
            [
                self._process_process_study_table_row(row)
                for row in self.study_table.collect()
            ],
        ).persist()

        # Adding study level metadata:
        evidence = (
            all_hits
            # Joining study level metadata:
            .join(self.study_table, on="studyId", how="inner")
            # Selecting relevant columns:
            .select(
                
                "diseaseFromSourceMappedId",
                "projectDescription",
                "projectId",
                "studyId",
                "studyOverview",
                "contrast",
                "crisprScreenLibrary",
                "cellType",
                "cellLineBackground",
                "geneticBackground",
                "statisticalTestTail",
                "resourceScore",
                "log2FoldChangeValue",
                "releaseDate",
                "releaseVersion",
                # Study level metadata:

                # Evidence level data:
                "targetFromSourceId",
                f.explode("diseaseFromSourceMappedId").alias("diseaseFromSourceMappedId"),
                # Static fields:
                f.lit("ot_crispr").alias("datasourceId"),
                f.lit("ot_partner").alias("datatypeId"),

            )

        ).
        return None


def main(study_table, output_file, data_folder) -> None:
    # Get logger:
    logger = logging.getLogger(__name__)

    # Initialize spark session:
    spark = initialize_sparksession()

    # Report input data:
    logger.info(f"Study information is read from: {study_table}")
    logger.info(f"Evidence saved to: {output_file}")
    logger.info(f"Data files read from: {data_folder}")

    # Parsing study table:
    study_table_raw = spark.read.csv(
        study_table, header=True, inferSchema=True, sep=","
    )
    study_table = otar_crispr_study_parser(study_table_raw)

    # # Parsing study table:
    # parser = OTAR_CRISPR_study_parser_old(study_table)

    # # Read data:
    # parser.generate_evidence(data_folder)

    # # Save data:
    # parser.write_evidence(output_file)


def parse_arguments() -> argparse.Namespace:
    """Parsing command line arguments.

    Returns:
        argparse.Namespace: _description_
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

    # Initialize logger:
    initialize_logger(__name__, args.log_file)

    study_table = args.study_table
    log_file = args.log_file
    data_folder = args.data_folder
    output_file = args.output

    # Configure logger:
    initialize_logger(__name__, log_file)

    # Process data:
    main(study_table, output_file, data_folder)
