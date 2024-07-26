#!/usr/bin/env python3
"""Parser for data submitted from the Validation Lab.

- The parser uses results from various screening experiments.
- Input files are defined in the ValidationLab_config.json configuration in the PPP-evidence-configuration repo
- Data files:
    - Cell line data with biomarker annotation
    - Result file with the measured values.
    - Sanger Model Passport data.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from functools import reduce
from typing import Callable

import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Column, DataFrame, SparkSession

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    read_ppp_config,
    write_evidence_strings,
)


class BiomarkerParser:
    """This class is responsible for parsing biomarker data.

    There are a handful of accepted biomarkers listed as dictionary keys. These values are expected to be
    column names in the biomarker input table. Then the column values are mapped to a standardised format.

    The output expected to collate all the biomarkers into a single column for each cell line. With the
    follwoing schema:

    root
    |-- cellName: string (nullable = true)
    |-- biomarkers: array (nullable = false)
    |    |-- element: struct (containsNull = false)
    |    |    |-- name: string (nullable = false)
    |    |    |-- description: string (nullable = false)
    """

    BIOMARKERMAPS = {
        "PAN": {
            "direct_mapping": {
                "CO": {"name": "PAN-CO", "description": "Pan-colorectal carcinoma"}
            }
        },
        "MS_status": {
            "direct_mapping": {
                "MSI": {"name": "MSI", "description": "Microsatellite unstable"},
                "MSS": {"name": "MSS", "description": "Microsatellite stable"},
            }
        },
        "PAM50_status": {
            "direct_mapping": {
                "Luminal B": {
                    "name": "Luminal B",
                    "description": "PAM50 status: luminal B",
                },
                "Luminal A": {
                    "name": "Luminal B",
                    "description": "PAM50 status: luminal A",
                },
                "Basal": {"name": "Basal", "description": "PAM50 status: basal"},
            }
        },
        "Hormone_status": {
            "direct_mapping": {
                "Hormone dependent": {
                    "name": "Hormone dependent",
                    "description": "Hormone dependency: hormone dependent",
                },
                "Hormone refractory (TNBC)": {
                    "name": "Hormone refractory (TNBC)",
                    "description": "Hormone dependency:  hormone refractory (triple-negative breast cancer)",
                },
            }
        },
        "CRIS_subtype": {
            "direct_mapping": {
                "A": {
                    "name": "CRIS-A",
                    "description": "mucinous, glycolytic, enriched for microsatellite instability or KRAS mutations.",
                },
                "B": {
                    "name": "CRIS-B",
                    "description": "TGF-β pathway activity, epithelial-mesenchymal transition, poor prognosis.",
                },
                "C": {
                    "name": "CRIS-C",
                    "description": "elevated EGFR signalling, sensitivity to EGFR inhibitors.",
                },
                "D": {
                    "name": "CRIS-D",
                    "description": "WNT activation, IGF2 gene overexpression and amplification.",
                },
                "E": {
                    "name": "CRIS-E",
                    "description": "Paneth cell-like phenotype, TP53 mutations.",
                },
                "?": {"name": "CRIS-?", "description": "CRIS subtype undetermined."},
            }
        },
        "KRAS_status": {
            "description": "KRAS mutation status: ",
            "name": "KRAS-",
        },
        "TP53_status": {
            "description": "TP53 mutation status: ",
            "name": "TP53-",
        },
        "APC_status": {
            "description": "APC mutation status: ",
            "name": "APC-",
        },
        "BRAF_status": {"description": "BRAF mutation status: ", "name": "BRAF-"},
    }

    @staticmethod
    def _get_biomarker_wrapper(biomarker_map: dict) -> Callable:
        @f.udf(
            t.StructType(
                [
                    t.StructField("name", t.StringType(), False),
                    t.StructField("description", t.StringType(), False),
                ]
            )
        )
        def wrapped_function(column_name: Column, biomarker: Column) -> dict | None:
            """This function returns a struct with the biomarker name and description."""
            # If the biomarker is not applied on the given cell type the value is zero:
            if biomarker == "0":
                return None

            # If the biomarker value has a direct mapping:
            if "direct_mapping" in biomarker_map[column_name]:
                try:
                    return biomarker_map[column_name]["direct_mapping"][biomarker]
                except KeyError:
                    logging.warning(
                        f"Could not find direct mapping for {column_name}:{biomarker}"
                    )
                    return None

            # If the value needs to be parsed:
            if biomarker == "wt":
                return {
                    "name": biomarker_map[column_name]["name"] + biomarker,
                    "description": biomarker_map[column_name]["description"]
                    + "wild type",
                }
            elif biomarker == "mut":
                return {
                    "name": biomarker_map[column_name]["name"] + biomarker,
                    "description": biomarker_map[column_name]["description"] + "mutant",
                }
            else:
                logging.warning(
                    f"Could not find direct mapping for {column_name}: {biomarker}"
                )
                return None

        return wrapped_function

    @classmethod
    def get_biomarkers(
        cls: type[BiomarkerParser], biomarker_df: DataFrame
    ) -> DataFrame:
        biomarkers_in_data = [
            biomarker
            for biomarker in cls.BIOMARKERMAPS.keys()
            if biomarker in biomarker_df.columns
        ]
        for biomarker in biomarkers_in_data:
            biomarker_df = biomarker_df.withColumn(
                biomarker,
                cls._get_biomarker_wrapper(cls.BIOMARKERMAPS)(
                    f.lit(biomarker), f.col(biomarker)
                ),
            )

        # The biomarker columns are unstacked into one single 'biomarkers' column:
        biomarker_unstack = f"""stack({len(biomarkers_in_data)}, {", ".join([f"'{x}', {x}" for x in biomarkers_in_data])}) as (biomarker_name, biomarker)"""

        return (
            biomarker_df
            # Selecting cell line name, cell line annotation and applyting the stacking expression:
            .select(
                f.col("Cell Line Name").alias("cellLineName"),
                f.expr(biomarker_unstack),
            )
            # Filter out all null biomarkers:
            .filter((f.col("biomarker").isNotNull()))
            # Grouping data by cell lines again:
            .groupBy("cellLineName")
            .agg(
                f.collect_list("biomarker").alias("biomarkerList"),
            )
            .persist()
        )


@dataclass
class ValidationLabEvidenceParser:
    projects: list
    assays: list
    sharedParemeters: dict

    def generate_evidence(
        self: ValidationLabEvidenceParser,
        spark: SparkSession,
        input_path: str,
        model_passport_file: str,
    ) -> DataFrame:
        # Get diease cell lines:
        disease_cell_lines = (
            spark.read.option("multiline", True)
            .csv(model_passport_file, header=True, sep=",", quote='"')
            .select(
                f.col("model_name").alias("cellLineName"),
                f.col("model_id").alias("id"),
                f.lower(f.col("tissue")).alias("tissueFromSource"),
            )
        )

        return (
            # Combining evidence from all validation projects
            reduce(
                lambda df1, df2: df1.unionByName(df2),
                [
                    # Generate evidence from a given project:
                    ValidationLabProjectParser(**project).parse_evidence(
                        spark, input_path, self.assays
                    )
                    for project in self.projects
                    if not project["excludeStudy"]
                ],
            )
            .join(disease_cell_lines, on="cellLineName", how="left")
            .select(
                # Adding shared columns:
                *[
                    f.lit(value).alias(colname)
                    for colname, value in self.sharedParemeters.items()
                ],
                # constants:
                "studyOverview",
                "releaseDate",
                "releaseVersion",
                # Evidence:
                "targetFromSourceId",
                "diseaseFromSourceMappedId",
                "diseaseFromSource",
                # Assessments:
                "resourceScore",
                "assessments",
                "assays",
                # Cell lines:
                f.array(
                    f.struct(
                        f.col("cellLineName").alias("name"),
                        f.col("id").alias("id"),
                        f.col("tissueFromSource").alias("tissue"),
                        f.col("tissueId"),
                    )
                ).alias("diseaseCellLines"),
                "biomarkerList",
                # Primary project:
                "primaryProjectHit",
                "primaryProjectId",
            )
        )


@dataclass
class ValidationLabProjectParser:
    diseaseFromSourceMappedId: str
    diseaseFromSource: str
    tissueId: str
    experimentDataFile: str
    biomarkerDataFile: str
    studyOverview: str
    primaryProjectId: str
    releaseVersion: str
    releaseDate: str
    excludeStudy: bool

    def parse_evidence(
        self: ValidationLabProjectParser,
        spark: SparkSession,
        input_path: str,
        assays: list,
    ) -> DataFrame:
        logger = logging.getLogger(__name__)

        validation_input_file = f"{input_path}/{self.experimentDataFile}"
        biomarker_file = f"{input_path}/{self.biomarkerDataFile}"

        logger.info(f"Reading project data: {validation_input_file}")
        logger.info(f"Reading biomarker data: {biomarker_file}")

        # Reading biomarker data:
        raw_biomarkers = spark.read.csv(biomarker_file, sep="\t", header=True)

        # Get parsed biomarkers:
        biomarkers = BiomarkerParser.get_biomarkers(raw_biomarkers)

        # Reading data:
        raw_experiment_data = spark.read.option("multiline", True).csv(
            validation_input_file, sep="\t", header=True
        )

        # First round of processing:
        processed_experiment = (
            raw_experiment_data
            # Read full evidence data:
            .select(
                # Extract target name:
                f.col("gene_ID").alias("targetFromSourceId"),
                # Extract cell-line name:
                f.col("cell_line_ID").alias("cellLineName"),
                # Extract VL assessments:
                f.split(f.col("OTVL_Assessment"), r"\n").alias(
                    "assessments"
                ),
                f.col("OTAR Primary Project Hit")
                .cast(t.BooleanType())
                .alias("primaryProjectHit"),
                # Extract assessment score:
                f.col("OTVL_Assessment_Score")
                .cast(t.FloatType())
                .alias("resourceScore"),
                # Extract all assays:
                *[
                    f.col(assay["label"])
                    .cast(t.BooleanType())
                    .alias(assay["shortName"])
                    for assay in assays
                    if assay["label"] in raw_experiment_data.columns
                ],
            )
        )

        # Formatting assays (depends on the formatted experiment data):
        parsed_assays = self._parse_assay(processed_experiment, assays, spark)

        # attributes added as column:
        columns_to_add = [
            "diseaseFromSourceMappedId",
            "diseaseFromSource",
            "tissueId",
            "studyOverview",
            "primaryProjectId",
            "releaseVersion",
            "releaseDate",
        ]

        return (
            processed_experiment.join(
                parsed_assays, on=["targetFromSourceId", "cellLineName"], how="inner"
            )
            .select(
                # Columns from the raw evidence file:
                "targetFromSourceId",
                "cellLineName",
                "assessments",
                "resourceScore",
                "primaryProjectHit",
                # Column from the assay parser:
                "assays",
                # Columns from the project metadata:
                *[
                    f.lit(getattr(self, colname)).alias(colname)
                    for colname in columns_to_add
                ],
            )
            .join(biomarkers, on="cellLineName", how="left")
        )

    @staticmethod
    def _parse_assay(
        raw_evidence_df: DataFrame, assays: list, spark: SparkSession
    ) -> DataFrame:
        """Organise experimental data into the right shape.

        For each evidence there might be a varying number of assays available.
        This function collate all the assays into a single list column of structs.
        """
        # Generate unpivot expression - not all assays are present in the data, some project might done with fewer assays:
        assay_names = [
            assay["shortName"]
            for assay in assays
            if assay["shortName"] in raw_evidence_df.columns
        ]

        unpivot_expression = f"""stack({len(assay_names)}, {', '.join([f"'{assay}', `{assay}`" for assay in assay_names])}) as (shortName, isHit)"""

        return (
            raw_evidence_df.select(
                "targetFromSourceId", "cellLineName", f.expr(unpivot_expression)
            )
            # Dropping rows without assessment:
            .filter(f.col("isHit").isNotNull())
            # Joining with the assay metadata:
            .join(spark.createDataFrame(assays), on="shortName", how="inner")
            .groupBy("targetFromSourceId", "cellLineName")
            .agg(
                f.collect_set(f.struct("shortName", "description", "isHit")).alias(
                    "assays"
                )
            )
        )


def main(
    config_file: str, output_file: str, cell_passport_file: str, data_folder: str
) -> None:
    # Logginging the configuration:
    logger = logging.getLogger(__name__)
    logger.info(f"Config file: {config_file}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"Sanger Model Passport file: {cell_passport_file}")
    logger.info(f"Validation lab files read from: {data_folder}")

    # Initialize spark session
    spark = initialize_sparksession()

    # Parse experimental parameters:
    validation_lab_config = read_ppp_config(config_file)
    logger.info(f"Number of projects: {len(validation_lab_config['projects'])}")
    logger.info(f"Number of assays: {len(validation_lab_config['assays'])}")
    logger.info(
        f"List of assays: {', '.join([assay['shortName'] for assay in validation_lab_config['assays']])}"
    )

    # Generate evidence from all validation projects:
    combined_evidence_df = ValidationLabEvidenceParser(
        **validation_lab_config
    ).generate_evidence(spark, data_folder, cell_passport_file)

    logger.info(
        f"Number of rows in the combined evidence: {combined_evidence_df.count()}"
    )
    # Save the combined evidence dataframe:
    write_evidence_strings(combined_evidence_df, output_file)
    logger.info(f"Saved evidence to {output_file}")


def parse_arguments() -> argparse.Namespace:
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="This script parse validation lab data and generates disease target evidence."
    )
    parser.add_argument(
        "--output_file", "-o", type=str, help="Output file. gzipped JSON", required=True
    )
    parser.add_argument(
        "--parameter_file",
        "-i",
        type=str,
        help="A JSON file describing exeriment metadata",
        required=True,
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="File into which the logs are saved",
        required=False,
    )
    parser.add_argument(
        "--cell_passport_file", type=str, help="Cell passport file", required=True
    )
    parser.add_argument(
        "--data_folder",
        type=str,
        help="Location of input files to process.",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    # Reading output file name from the command line:
    args = parse_arguments()
    initialize_logger(args.log_file)

    # Passing all the required arguments:
    main(
        args.parameter_file, args.output_file, args.cell_passport_file, args.data_folder
    )
