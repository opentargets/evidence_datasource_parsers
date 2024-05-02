#!/usr/bin/env python3
"""Parser for data submitted from the Validation Lab.

- The parser uses results from various screening experiments.
- Input files are defined in the ValidationLab_config.json configuration in the PPP-evidence-configuration repo
- Files:
    - Cell line data with biomarker annotation
    - Result file with the measured values.
    - Hypothesis files with the binary flags for each tested biomarker in the context of the tested gene.
"""

import argparse
import logging
import sys
from dataclasses import dataclass
from functools import reduce

import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    read_ppp_config,
    write_evidence_strings,
)


class BiomarkerParser:
    BIOMARKERMAPS = {
        "PAN": {
            "direct_mapping": {
                "CO": {"name": "PAN-CO", "description": "Pan-colorectal carcinoma"}
            }
        },
        "MS_status": {
            "direct_mapping": {
                "MSI": {"name": "MSI", "description": "Microsatellite instable"},
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
    @f.udf(
        t.StructType(
            [
                t.StructField("name", t.StringType(), False),
                t.StructField("description", t.StringType(), False),
            ]
        )
    )
    def _get_biomarker(column_name, biomarker, BIOMARKERMAPS):
        """This function returns a struct with the biomarker name and description."""

        # If the biomarker has a direct mapping:
        if "direct_mapping" in BIOMARKERMAPS[column_name]:
            try:
                return BIOMARKERMAPS[column_name]["direct_mapping"][biomarker]
            except KeyError:
                logging.warning(
                    f"Could not find direct mapping for {column_name}:{biomarker}"
                )
                return None

        # If the value needs to be parsed:
        if biomarker == "wt":
            return {
                "name": BIOMARKERMAPS[column_name]["name"] + biomarker,
                "description": BIOMARKERMAPS[column_name]["description"] + "wild type",
            }
        elif biomarker == "mut":
            return {
                "name": BIOMARKERMAPS[column_name]["name"] + biomarker,
                "description": BIOMARKERMAPS[column_name]["description"] + "mutant",
            }
        else:
            logging.warning(
                f"Could not find direct mapping for {column_name}: {biomarker}"
            )
            return None

    @classmethod
    def get_biomarkers(cls, biomarker_df):
        biomarkers_in_data = [
            biomarker
            for biomarker in cls.BIOMARKERMAPS.keys()
            if biomarker in biomarker_df.columns
        ]

        # Applying the full map on the dataframe one-by-one:
        processed_biomarkers = reduce(
            lambda DF, value: DF.withColumn(*value),
            map(
                # Function to process biomarker:
                lambda biomarker: (
                    biomarker,
                    cls._get_biomarker(
                        f.lit(biomarker), f.col(biomarker), cls.BIOMARKERMAPS
                    ),
                ),
                # Iterator to apply the function over:
                biomarkers_in_data,
            ),
            biomarker_df,
        )

        # The biomarker columns are unstacked into one single 'biomarkers' column:
        biomarker_unstack = f"""stack({len(biomarkers_in_data)}, {", ".join([f"'{x}', {x}" for x in biomarkers_in_data])}) as (biomarker_name, biomarkers)"""

        return (
            processed_biomarkers
            # Selecting cell line name, cell line annotation and applyting the stacking expression:
            .select(
                f.col("Cell Line Name").alias("cellLineName"),
                f.expr(biomarker_unstack),
            )
            # Filter out all null biomarkers:
            .filter((f.col("biomarkers").isNotNull()))
            # Grouping data by cell lines again:
            .groupBy("cellLineName")
            .agg(
                f.collect_list("biomarkers").alias("biomarkers"),
            )
            .persist()
        )


@dataclass
class ValidationLabEvidenceParser:
    projects: list
    assays: list
    sharedParemeters: dict

    @staticmethod
    def _get_disease_cell_lines(model_passport_df: DataFrame) -> DataFrame:
        return (
            # The following option is required to correctly parse CSV records which contain newline characters.
            model_passport_df.select(
                f.col("model_name").alias("cellLineName"),
                f.col("model_id").alias("id"),
                f.lower(f.col("tissue")).alias("tissueFromSource"),
            )
        )

    def generate_evidence(
        self, spark: SparkSession, input_path: str, model_passport_file: str
    ):
        # Get diease cell lines:
        disease_cell_lines = self._get_disease_cell_lines(
            spark.read.option("multiline", True).csv(
                model_passport_file, header=True, sep=",", quote='"'
            )
        )

        return (
            # Combining evidence from all validation projects
            reduce(
                lambda df1, df2: df1.unionByName(df1),
                [
                    # Generate evidence from a given project:
                    ValidationLabProjectParser(**project).parse_evidence(
                        spark, input_path, self.assays
                    )
                    for project in self.projects
                    if not project["excludeStudy"]
                ],
            )
            # Adding shared columns:
            .select(
                "*",
                *[
                    f.lit(value).alias(colname)
                    for colname, value in self.sharedParemeters.items()
                ],
            )
            .join(disease_cell_lines, on="cellLineName", how="left")
            .select(
                # constants:
                "datasourceId",
                "datatypeId",
                "studyOverview",
                "projectId",
                "releaseDate",
                "releaseVersion",
                # Evidence:
                "targetFromSource",
                "diseaseFromSourceMappedId",
                "diseaseFromSource",
                # Assessments:
                f.lit(1.0).alias("resourceScore"),  # Resource score will come from VL.
                "assessment",
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
                "biomarkers",
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
        self, spark: SparkSession, input_path: str, assays: list
    ) -> DataFrame:
        # Reading data:
        raw_experiment_data = spark.read.option("multiline", True).csv(
            f"{input_path}/{self.experimentDataFile}", sep="\t", header=True
        )
        # First round of processing:
        processed_experiment = self._process_raw_evidence(raw_experiment_data, assays)

        # Reading biomarker data:
        raw_biomarkers = spark.read.csv(
            f"{input_path}/{self.biomarkerDataFile}", sep="\t", header=True
        )

        # Formatting assays (depends on the formatted experiment data):
        parsed_assays = self._parser_assay_object(processed_experiment, assays)

        # Get parsed biomarkers:
        biomarkers = BiomarkerParser.get_biomarkers(raw_biomarkers)

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
                parsed_assays, on=["targetFromSource", "cellLineName"], how="inner"
            )
            .select(
                # Columns from the raw evidence file:
                "targetFromSource",
                "cellLineName",
                "assessment",
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
    def _process_raw_evidence(raw_data: DataFrame, assays: list) -> DataFrame:
        # Read full evidence data:
        return raw_data.select(
            # Extract target name:
            f.col("gene_ID").alias("targetFromSource"),
            # Extract cell-line name:
            f.col("cell_line_ID").alias("cellLineName"),
            # Extract VL assessment:
            f.regexp_replace(f.col("OTVL Assessment"), r"\n", " ").alias("assessment"),
            f.col("OTAR Primary Project Hit")
            .cast(t.BooleanType())
            .alias("primaryProjectHit"),
            # Extract all assays:
            *[
                f.col(assay["label"]).cast(t.BooleanType()).alias(assay["shortName"])
                for assay in assays
            ],
        )

    @staticmethod
    def _parser_assay_object(raw_evidence_df: DataFrame, assays: list) -> DataFrame:
        """Organise experimental data into the right shape."""

        # Generate unpivot expression:
        assay_names = [assay["shortName"] for assay in assays]

        unpivot_expression = f"""stack({len(assay_names)}, {', '.join([f"'{assay}', `{assay}`" for assay in assay_names])}) as (shortName, isAssayHit)"""

        return (
            raw_evidence_df.select(
                "targetFromSource", "cellLineName", f.expr(unpivot_expression)
            )
            .join(spark.createDataFrame(assays), on="shortName", how="inner")
            .groupBy("targetFromSource", "cellLineName")
            .agg(f.collect_set(f.struct("shortName", "description")).alias("assays"))
        )


def main(
    config_file: str, output_file: str, cell_passport_file: str, data_folder: str
) -> None:
    # Initialize spark session
    spark = initialize_sparksession()

    # Parse experimental parameters:
    parameters = read_ppp_config(config_file)

    # Opening and parsing the cell passport data from Sanger:
    cell_passport_df = get_cell_passport_data(spark, cell_passport_file)

    logging.info(f"Cell passport dataframe has {cell_passport_df.count()} rows.")

    logging.info("Parsing experiment data...")

    # Create evidence for all experiments:
    evidence_dfs = []
    for experiment in parameters["experiments"]:
        evidence_dfs.append(
            parse_experiment(spark, experiment, cell_passport_df, data_folder)
        )

    # combine all evidence dataframes into one:
    combined_evidence_df = reduce(lambda df1, df2: df1.union(df2), evidence_dfs)

    # Save the combined evidence dataframe:
    write_evidence_strings(combined_evidence_df, output_file)


if __name__ == "__main__":
    # Reading output file name from the command line:
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
    args = parser.parse_args()

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Passing all the required arguments:
    main(
        args.parameter_file, args.output_file, args.cell_passport_file, args.data_folder
    )
