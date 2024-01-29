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
from functools import reduce

from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql.dataframe import DataFrame

from common.evidence import (
    write_evidence_strings,
    initialize_sparksession,
    read_ppp_config,
)


# Datasource-wide constants:
VALIDATION_LAB_DATASOURCE_ID = "ot_crispr_validation"
VALIDATION_LAB_DATATYPE_ID = "ot_validation_lab"

# This is a map that provides recipie to generate the biomarker objects
# If a value cannot be found in the map, the value will be returned.
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


class ParseHypotheses:
    def __init__(self, spark) -> None:
        self.spark = spark

    def parse_hypotheses(self, expectedFile: str, observedFile: str) -> DataFrame:
        """
        Hypothesis is parsed from two files describing the expected and observed results.
        This function reads the files, compare them, parses the hypothesis as biomarker? + status

        Args:
            expectedFile: file with the expected results
            observedFile: file with the observed results
        Returns:
            DataFrame with the following schema:

            |-- gene: string (nullable = true)
            |-- hypotheses: array (nullable = false)
            |    |-- element: struct (containsNull = false)
            |    |    |-- name: string (nullable = true)
            |    |    |-- description: string (nullable = true)
            |    |    |-- status: string (nullable = false)
        """

        # The observed and expected results follows the same schema and parsed the same way:
        expected_df = self.read_hypothesis_data(expectedFile, "expected")
        observed_df = self.read_hypothesis_data(observedFile, "observed")

        return (
            expected_df
            # Joining expected vs observed hypothesis tables:
            .join(observed_df, on=["gene", "hypothesis"], how="inner")
            # Filter hypotheses where at least one was True:
            # .filter(col('expected') | col('observed'))
            # From the hypothesis column eg. CRIS_subtype-B ectract the type CRIS_subtype and the call: B
            .withColumn(
                "hypothesis_type", f.element_at(f.split(f.col("hypothesis"), "-"), 1)
            )
            .withColumn(
                "hypothesis_call", f.element_at(f.split(f.col("hypothesis"), "-"), 2)
            )
            # Using the biomarker parser generate an struct similar to the biomarker object:
            .withColumn(
                "hypothesis",
                get_biomarker(f.col("hypothesis_type"), f.col("hypothesis_call")),
            )
            # Besides the annotation we add the status of the hypothesis:
            .withColumn(
                "status",
                f.when(f.col("expected") & f.col("observed"), "observed and expected")
                .when(f.col("expected"), "expected but not observed")
                .when(f.col("observed"), "observed but not expected")
                .otherwise("not expected and not observed"),
            )
            .withColumn("hypothesis", f.struct("hypothesis.*", "status"))
            # Collect all hypotheses for each gene:
            .groupBy("gene")
            .agg(f.collect_set("hypothesis").alias("validationHypotheses"))
            .persist()
        )

    def read_hypothesis_data(self, file: str, call: str) -> DataFrame:
        """Parsing the hypothesis file.

        Args:
            file: hypothesis file tsv with genes in rows and biomarkers in columns, the hypotheses are boolean.
        Returns:
            DataFrame with the following columns: gene, hypothesis, call (true/false)
        """

        hypothesis_df = (
            self.spark.read.csv(file, sep="\t", header=True).withColumnRenamed(
                "Gene", "gene"
            )
            # The gene names are manually typed, there are lower and upper case names:
            .withColumn("gene", f.upper(f.col("gene")))
        )

        # The first column is the gene name, the rest are the hypotheses:
        hypothesis_columns = hypothesis_df.columns[1:]

        unpivot_expression = f"""stack({len(hypothesis_columns)}, {", ".join([f"'{x}', `{x}`" for x in hypothesis_columns])} ) as (hypothesis, {call})"""

        return (
            hypothesis_df.select("Gene", f.expr(unpivot_expression))
            .withColumn(call, f.col(call).cast(t.BooleanType()))
            .persist()
        )


@f.udf(
    t.StructType(
        [
            t.StructField("name", t.StringType(), False),
            t.StructField("description", t.StringType(), False),
        ]
    )
)
def get_biomarker(column_name, biomarker):
    """This function returns with a struct with the biomarker name and description"""

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
        logging.warning(f"Could not find direct mapping for {column_name}: {biomarker}")
        return None


def get_cell_passport_data(spark: SparkSession, cell_passport_file: str) -> DataFrame:

    # loading cell line annotation data from Sanger:
    return (
        spark.read.option("multiline", True)
        .csv(cell_passport_file, header=True, sep=",", quote='"')
        .select(
            f.regexp_replace(f.col("model_name"), r"-", "").alias("cellName"),
            f.col("model_id").alias("cellId"),
            f.col("tissue"),
        )
        .persist()
    )


def parse_experiment(
    spark: SparkSession, parameters: dict, cellPassportDf: DataFrame, data_folder: str
) -> DataFrame:
    """
    Parse experiment data from a file.

    Args:
        spark: Spark session.
        parameters: Dictionary of experimental parameters.
        cellPassportDf: Dataframe of cell passport data.
        data_folder: Location of the input data files.

    Returns:
        A dataframe of experiment data.
    """

    # Extracting parameters:
    experiment_file = f"{data_folder}/{parameters['experimentData']}"
    contrast = parameters["contrast"]
    studyOverview = parameters["studyOverview"]
    projectId = parameters["projectId"]
    projectDescription = parameters["projectDescription"]
    diseaseFromSource = parameters["diseaseFromSource"]
    diseaseFromSourceMapId = parameters["diseaseFromSourceMappedId"]
    confidenceCutoff = parameters["confidenceCutoff"]
    cell_line_file = f"{data_folder}/{parameters['cellLineFile']}"
    tissue_id = parameters["tissueId"]

    # The hypothesis is defined by two datasets:
    hypothesis_expected_file = f"{data_folder}/{parameters['hypothesisFileExpected']}"
    hypothesis_observed_file = f"{data_folder}/{parameters['hypothesisFileObserved']}"

    # Reading cell metadata from validation lab:
    validation_lab_cell_lines = (
        spark.read.csv(cell_line_file, sep="\t", header=True)
        # Renaming columns:
        .withColumnRenamed("cell_line", "cellName")
        .drop("tissue")
        # Joining dataset with cell model data read downloaded from Sanger website:
        .join(cellPassportDf, on="cellName", how="left")
        # Adding UBERON code to tissues (it's constant colon)
        .withColumn("tissueID", f.lit(tissue_id))
        # generating disease cell lines object:
        .withColumn(
            "diseaseCellLines",
            f.array(
                f.struct(
                    f.col("cellName").alias("name"),
                    f.col("cellId").alias("id"),
                    f.col("tissue"),
                    f.lit(tissue_id).alias("tissueId"),
                )
            ),
        )
        .persist()
    )

    logging.info(
        f"Validation lab cell lines has {validation_lab_cell_lines.count()} cell types."
    )

    # Defining how to process biomarkers:
    # 1. Looping through all possible biomarker - from biomarkerMaps.keys()
    # 2. The biomakers are then looked up in the map and process based on how the map defines.
    # 3. Description is also added read from the map.
    biomarkers_in_data = [
        biomarker
        for biomarker in BIOMARKERMAPS.keys()
        if biomarker in validation_lab_cell_lines.columns
    ]

    expressions = (
        # Function to process biomarker:
        lambda biomarker: (
            biomarker,
            get_biomarker(f.lit(biomarker), f.col(biomarker)),
        ),
        # Iterator to apply the function over:
        biomarkers_in_data,
    )

    # Applying the full map on the dataframe one-by-one:
    biomarkers = reduce(
        lambda DF, value: DF.withColumn(*value), expressions, validation_lab_cell_lines
    )

    # The biomarker columns are unstacked into one single 'biomarkers' column:
    biomarker_unstack = f"""stack({len(biomarkers_in_data)}, {", ".join([f"'{x}', {x}" for x in biomarkers_in_data])}) as (biomarker_name, biomarkers)"""

    validation_lab_cell_lines = (
        biomarkers
        # Selecting cell line name, cell line annotation and applyting the stacking expression:
        .select(
            f.col("cellName"),
            "diseaseCellLines",
            f.expr(biomarker_unstack),
        )
        # Filter out all null biomarkers:
        .filter(
            (f.col("biomarkers").isNotNull())
            &
            # Following the request of the validation lab, we are removing CRIS biomarker annotation:
            (f.col("biomarker_name") != "CRIS_subtype")
        )
        # Grouping data by cell lines again:
        .groupBy("cellName")
        .agg(
            f.collect_list("biomarkers").alias("biomarkers"),
            f.first(f.col("diseaseCellLines")).alias("diseaseCellLines"),
        )
        .persist()
    )

    # Reading and processing hypothesis data:
    hypothesis_generator = ParseHypotheses(spark)
    validation_hypotheses_df = hypothesis_generator.parse_hypotheses(
        expectedFile=hypothesis_expected_file, observedFile=hypothesis_observed_file
    )

    # Reading experiment data from validation lab:
    evidence = (
        # Reading evidence:
        spark.read.csv(experiment_file, sep="\t", header=True)
        .withColumnRenamed("cell_line", "cellName")
        # Genes need to be uppercase:
        .withColumn("gene", f.upper(f.col("gene")))
        # Joining hypothesis data:
        .join(validation_hypotheses_df, on="gene", how="left")
        # Joining with cell line data:
        .join(validation_lab_cell_lines, on="cellName", how="left")
        # Selecting all columns:
        .select(
            f.col("gene").alias("targetFromSourceId"),
            f.col("validationHypotheses"),
            f.when(
                f.col("effect_size").cast("double") > 0,
                f.col("effect_size").cast("double"),
            )
            .otherwise(0)
            .alias("resourceScore"),
            f.when(
                f.col("effect_size").cast("double") >= confidenceCutoff,
                f.lit("significant"),
            )
            .otherwise(f.lit("not significant"))
            .alias("confidence"),
            f.when(f.col("expected_to_pass") == "TRUE", f.lit("significant"))
            .otherwise(f.lit("not significant"))
            .alias("expectedConfidence"),
            f.lit("upper tail").alias("statisticalTestTail"),
            f.lit(contrast).alias("contrast"),
            f.lit(studyOverview).alias("studyOverview"),
            f.lit(diseaseFromSourceMapId).alias("diseaseFromSourceMappedId"),
            f.lit(diseaseFromSource).alias("diseaseFromSource"),
            f.lit(projectId).alias("projectId"),
            f.lit(projectDescription).alias("projectDescription"),
            f.col("biomarkers").alias("biomarkerList"),
            f.col("diseaseCellLines"),
            f.lit(VALIDATION_LAB_DATATYPE_ID).alias("datatypeId"),
            f.lit(VALIDATION_LAB_DATASOURCE_ID).alias("datasourceId"),
        )
        .persist()
    )

    logging.info(f"Evidence count: {evidence.count()}.")
    return evidence


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
