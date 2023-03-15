#!/usr/bin/env python3
"""Parser for data submitted by the Encore team."""

import argparse
import json
import logging
import os
import sys
from functools import reduce

from scipy import stats as st

from typing import Dict, List, Any

from pyspark.sql import SparkSession, functions as f, types as t
from pyspark.sql.dataframe import DataFrame
from pyspark.sql.window import Window

from common.evidence import (
    initialize_sparksession,
    write_evidence_strings,
    GenerateDiseaseCellLines,
)


class EncoreEvidenceGenerator:
    """This parser generates disease/target evidence based on datafiles submitted by the Encore team.

    Input for the parser:
        - JSON configuration describing the experiments.

    For each experiment the following input files are expected:
        - lfc_file: log fold change values for every gene pairs, for every experiment.
        - bliss_file: the analysed cooperative effect for each gene pairs calculated by the BLISS analysis.
        - gemini_file: the analysed cooperative effect for each gene pairs calculated by the Gemini analysis.

    The parser joins the lfc data with the cooperativity measurements by gene pairs and cell line.
    Apply filter on gene pairs:
        - only extract gene pairs with significant log fold change.
        - only extract gene pairs with significant cooperative effect.
    """

    def __init__(
        self, spark: SparkSession, cell_passport_df: DataFrame, shared_parameters: dict
    ) -> None:
        self.spark = spark
        self.cell_passport_df = cell_passport_df

        # Parsing paramters:
        self.log_fold_change_cutoff_p_val = shared_parameters["logFoldChangeCutoffPVal"]
        self.log_fold_change_cutoff_FDR = shared_parameters["logFoldChangeCutoffFDR"]
        self.interaction_cutoff_p_val = shared_parameters["blissCutoffPVal"]
        self.data_folder = shared_parameters["data_folder"]
        self.release_version = shared_parameters["releaseVersion"]
        self.release_date = shared_parameters["releaseDate"]

    @staticmethod
    @f.udf(
        t.ArrayType(
            t.StructType(
                [
                    t.StructField("targetFromSourceId", t.StringType(), False),
                    t.StructField("targetRole", t.StringType(), False),
                    t.StructField(
                        "interactingTargetFromSourceId", t.StringType(), False
                    ),
                    t.StructField("interactingTargetRole", t.StringType(), False),
                ]
            )
        )
    )
    def parse_targets(gene_pair: str, gene_role: str) -> List[Dict[Any, Any]]:
        """The gene pair string is split and assigned to the relevant role + exploding into evidence of two targets.

        gene pair: 'SHC1~ADAD1', where 'SHC1' is the library gene, and 'ADAD1' is the anchor gene.

        Both genes will be targetFromSource AND interactingTargetFromSource, while keeping their roles.
        """
        genes = gene_pair.split("~")
        roles = [gene_role.replace("Combinations", "").lower(), "anchor"]

        assert len(genes) == 2
        parsed = []

        for i, (gene, role) in enumerate(zip(genes, roles)):
            parsed.append(
                {
                    "targetFromSourceId": gene,
                    "targetRole": role,
                    "interactingTargetFromSourceId": genes[1] if i == 0 else genes[0],
                    "interactingTargetRole": roles[1] if i == 0 else roles[0],
                }
            )

        return parsed

    def get_lfc_data(self, lfc_file: str) -> DataFrame:
        """
        The table has results from CRISPR/Cas9 experiments: the p-value, false discovery rate (fdr)
        and the log fold-change file describes the observed survival of the cells.

        Each row is a gene pair, columns contain results from the measurements. Column names contain
        information on the tested cell line (SIDM) and ocassionally the replicate (CSID) identifiers.

        This really wide table is stacked to get only three numeric columns (p-value, fold-change and FDR),
        plus string column with cell line.

        Args:
          lfc_file (str): str

        Returns:
          A dataframe with the following columns:
            id
            cellLineName
            Note1
            Note2
            phenotypicConsequenceLogFoldChange
            phenotypicConsequencePValue
            phenotypicConsequenceFDR
        """

        # Fixed statistical field names:
        stats_fields = ["p-value", "fdr", "lfc"]

        # Reading the data into a single dataframe:
        lfc_df = self.spark.read.csv(lfc_file, sep="\t", header=True)

        # Collect the cell lines from the lfc file header:
        cell_lines = set(["_".join(x.split("_")[:-1]) for x in lfc_df.columns[4:]])

        # Generating struct for each cell lines:
        # SIDM are Sanger model identifiers, while CSID are replicate identifiers.
        # SIDM00049_CSID1053_p-value
        # SIDM00049_CSID1053_fdr
        # SIDM00049_CSID1053_lfc
        # Into: SIDM00049_CSID1053: struct(p-value, fdr, lfc)
        expressions = map(
            lambda cell: (
                cell,
                f.struct(
                    [
                        f.col(f"{cell}_{x}").cast(t.FloatType()).alias(x)
                        for x in stats_fields
                    ]
                ),
            ),
            cell_lines,
        )

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, lfc_df)

        # Stack the previously generated columns:
        unpivot_expression = f"""stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)"""

        return (
            res_df
            # Unpivot:
            .select("id", "Note1", "Note2", f.expr(unpivot_expression))
            # Extracting and renaming log-fold-change parameters:
            .select(
                "id",
                f.split(f.col("cellLineName"), "_").getItem(0).alias("cellLineName"),
                "Note1",
                "Note2",
                f.col("cellLineData.lfc").alias("phenotypicConsequenceLogFoldChange"),
                f.col("cellLineData.p-value").alias("phenotypicConsequencePValue"),
                f.col("cellLineData.fdr").alias("phenotypicConsequenceFDR"),
            )
        )

    def get_bliss_data(self, bliss_file: str) -> DataFrame:
        """
        It reads a bliss file, extracts the cell lines from the header, generates a struct for each cell
        line, unpivots the data, and finally extracts the bliss score and p-value

        Args:
          bliss_file (str): The path to the bliss file.

        Returns:
          A dataframe with the following columns:
            - id
            - cellLineName
            - geneticInteractionScore
            - geneticInteractionPValue
            - statisticalMethod
        """

        stats_fields = [
            "gene1",
            "gene2",
            "observed",
            "observed_expected",
            "pval",
            "zscore",
        ]
        # Read bliss file:
        bliss_df = self.spark.read.csv(bliss_file, sep="\t", header=True)

        # Collect cell-line/recplicate pairs from the headers:
        cell_lines = set(
            [
                "_".join(x.split("_")[0:2])
                for x in bliss_df.columns[4:]
                if x.startswith("SID")
            ]
        )

        # Generating struct for each cell lines:
        expressions = map(
            lambda cell: (
                cell,
                f.struct([f.col(f"{cell}_{x}").alias(x) for x in stats_fields]),
            ),
            cell_lines,
        )

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, bliss_df)

        # Stack the previously generated columns:
        unpivot_expression = f"""stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)"""

        # Windowing over the cell-line/gene-pair pairs:
        w = Window.partitionBy(["id", "cellLineName"]).orderBy(
            f.asc("geneticInteractionPValue")
        )

        return (
            res_df.select(
                # Create a consistent id column:
                f.regexp_replace(f.col("Gene_Pair"), ";", "~").alias("id"),
                # Unpivot:
                f.expr(unpivot_expression),
            )
            # Extracting and renaming bliss statistical values:
            .select(
                "id",
                f.split(f.col("cellLineName"), "_")[0].alias("cellLineName"),
                f.split(f.col("cellLineName"), "_")[1].alias("experimentId"),
                f.col("cellLineData.zscore").cast(t.FloatType()).alias("zscore"),
            )
            .filter(f.col("zscore") != "NaN")
            # BLISS data is not aggregated for cell lines, instead we have data for each replicates.
            # To allow averaging, data needs to be grouped by gene pair and cell line:
            .groupby("id", "cellLineName")
            # Averaging z-scores:
            .agg(
                f.sum(f.col("zscore")).alias("sumzscore"),
                f.count(f.col("zscore")).alias("experiment_count"),
            )
            .withColumn(
                "geneticInteractionScore",
                f.col("sumzscore") / f.sqrt(f.col("experiment_count")),
            )
            # Calculating p-value for the averaged z-score:
            .withColumn(
                "geneticInteractionPValue",
                f.udf(
                    lambda zscore: float(st.norm.sf(abs(zscore)) * 2)
                    if zscore
                    else None,
                    t.FloatType(),
                )(f.col("geneticInteractionScore")),
            )
            .select(
                "id",
                "cellLineName",
                f.lit("bliss").alias("statisticalMethod"),
                "geneticInteractionPValue",
                "geneticInteractionScore",
                # Based on the z-score we can tell the nature of the interaction:
                # - positive z-score means synergictic
                # - negative z-score means antagonistics interaction
                f.when(f.col("geneticInteractionScore") > 0, "synergistic")
                .otherwise("antagonistic")
                .alias("geneInteractionType"),
            )
        )

    def get_gemini_data(self, gemini_file: str) -> DataFrame:
        """Parsing cooperative effects calculated by gemini method

        The structure of the file similar to the lfc file, the process follows a similar fashion,
        except it slightly different. Different enough to get a separate function.
        """

        # Fixed statistical field names:
        stats_fields = ["score", "pval", "FDR"]

        # Reading the data into a single dataframe:
        gemini_df = self.spark.read.csv(gemini_file, sep="\t", header=True)

        # Collect the cell lines from the lfc file header:
        cell_lines = set(
            [
                "_".join(x.split("_")[:-1])
                for x in gemini_df.columns[4:]
                if x.startswith("SID")
            ]
        )

        # There are some problems in joining gemini files on Encore side. It causes a serious issues:
        # 1. Multiple Gene_Pair columns in the file -> these will be indexed in the pyspark dataframe
        # 2. Some columns for some cell lines will be missing eg. pvalue for SIDM00049_CSID1053
        #
        # To mitigate these issue we have to check for gene pair header and remove cell lines with incomplete data.
        if "Gene_Pair" in gemini_df.columns:
            gene_column = "Gene_Pair"
        elif "Gene_Pair0" in gemini_df.columns:
            gene_column = "Gene_Pair0"
        else:
            raise ValueError(
                f'No \'Gene_Pair\' column in Gemini data: {",".join(gemini_df.columns)}'
            )

        # We check if all stats columns available for all cell lines (this coming from a data joing bug at encore):
        missing_columns = [
            f"{cell}_{stat}"
            for cell in cell_lines
            for stat in stats_fields
            if f"{cell}_{stat}" not in gemini_df.columns
        ]
        cells_to_drop = set(["_".join(x.split("_")[:-1]) for x in missing_columns])

        # If there are missingness, the relevant cell lines needs to be removed from the analysis:
        if missing_columns:
            logging.warn(f'Missing columns: {", ".join(missing_columns)}')
            logging.warn(f'Dropping cell_lines: {", ".join(cells_to_drop)}')

            # Removing missing cell lines:
            cell_lines = [x for x in cell_lines if x not in cells_to_drop]

        # Generating struct for each cell lines:
        expressions = map(
            lambda cell: (
                cell,
                f.struct([f.col(f"{cell}_{x}").alias(x) for x in stats_fields]),
            ),
            cell_lines,
        )

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, gemini_df)

        # Stack the previously generated columns:
        unpivot_expression = f"""stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)"""

        return (
            res_df
            # Create a consistent id column:
            .withColumn("id", f.regexp_replace(f.col(gene_column), ";", "~"))
            # Unpivot:
            .select("id", f.expr(unpivot_expression))
            # Extracting and renaming gemini statistical values:
            .select(
                "id",
                f.regexp_replace("cellLineName", "_strong", "").alias("cellLineName"),
                f.col("cellLineData.score")
                .cast(t.FloatType())
                .alias("geneticInteractionScore"),
                f.col("cellLineData.pval")
                .cast(t.FloatType())
                .alias("geneticInteractionPValue"),
                f.col("cellLineData.FDR")
                .cast(t.FloatType())
                .alias("geneticInteractionFDR"),
                f.lit("gemini").alias("statisticalMethod"),
            )
        )

    def parse_experiment(self, parameters: dict) -> DataFrame:
        """Parsing experiments based on the experimental descriptions.

        Args:
            parameters: Dictionary of experimental parameters. The following keys are required:
                - dataset: Name of the dataset eg. COLO1- referring to the first libraryset of colo.
                - diseaseFromSource: Name of the disease model of the experiment.
                - diseaseFromSourceMappedId: EFO ID of the disease model of the experiment.
                - logFoldChangeFile: File path to the log fold change file.
                - geminiFile: File path to the gemini file.
                - blissFile: File path to the bliss file.

        Returns:
            A pyspark dataframe of experiment data.

        Process:
            - Reading all files.
            - Joining files.
            - Applying filters.
                - Apply filter on lf change
                - Apply filter on cooperativity
                - Apply filter on target role (keeping only library genes, controls are dropped)
            - Adding additional columns + finalizing evidence model
        """

        disease_from_source = parameters["diseaseFromSource"]
        disease_from_source_mapped_id = parameters["diseaseFromSourceMappedId"]
        dataset = parameters["dataset"]
        log_fold_change_file = parameters["logFoldChangeFile"]
        # gemini_file = parameters["geminiFile"] # <= Removed as genini method is dropped for now.
        blissFile = parameters["blissFile"]

        # Testing if experiment needs to be skipped:
        if parameters["skipStudy"] is True:
            logging.info(f"Skipping study: {dataset}")
            return None

        logging.info(f"Parsing experiment: {dataset}")

        # if no log fold change file is provided, we will not generate any evidence.
        if log_fold_change_file is None:
            logging.warning(f"No log fold change file provided for {dataset}.")
            return None

        # Reading lfc data:
        lfc_file = f"{self.data_folder}/{log_fold_change_file}"
        lfc_df = self.get_lfc_data(lfc_file)

        logging.info(
            f'Number of gene pairs in the log(fold change) dataset: {lfc_df.select("id").distinct().count()}'
        )
        logging.info(
            f'Number cell lines in the log(fold change) dataset: {lfc_df.select("cellLineName").distinct().count()}'
        )
        # A decistion was made to exclude gemini asessments, for now.
        # # Reading gemini data:
        # gemini_file = f"{self.data_folder}/{gemini_file}"
        # gemini_df = self.get_gemini_data(gemini_file)
        # logging.info(
        #     f'Number of gene pairs in the gemini dataset: {gemini_df.select("id").distinct().count()}'
        # )
        # logging.info(
        #     f'Number cell lines in the gemini dataset: {gemini_df.select("cellLineName").distinct().count()}'
        # )

        bliss_file = f"{self.data_folder}/{blissFile}"
        bliss_df = self.get_bliss_data(bliss_file)

        # Merging lfc + gemini:
        merged_dataset = (
            lfc_df
            # Data is joined by the gene-pair and cell line:
            .join(bliss_df, how="inner", on=["id", "cellLineName"])
            # Applying filters on logFoldChange + interaction p-value thresholds:
            .filter(
                (
                    f.col("phenotypicConsequencePValue")
                    <= self.log_fold_change_cutoff_p_val
                )
                & (f.col("phenotypicConsequenceFDR") <= self.log_fold_change_cutoff_FDR)
                & (f.col("geneticInteractionPValue") <= self.interaction_cutoff_p_val)
            )
            # Cleaning the cell line annotation:
            .withColumn("cellId", f.split(f.col("cellLineName"), "_").getItem(0))
            # Joining with cell passport data containing diseaseCellLine and biomarkers info:
            .join(
                self.cell_passport_df.select(
                    f.col("id").alias("cellId"), "diseaseCellLine", "biomarkerList"
                ),
                on="cellId",
                how="left",
            )
        )
        logging.info(
            f'Number of gene pairs in the merged dataset: {merged_dataset.select("id").count()}'
        )
        logging.info(
            f'Number of cell lines in the merged dataset: {merged_dataset.select("cellLineName").distinct().count()}'
        )

        evidence_df = (
            merged_dataset
            # Parsing and exploding gene names and target roles:
            .withColumn("id", self.parse_targets(f.col("id"), f.col("Note1")))
            .select("*", f.explode(f.col("id")).alias("genes"))
            .select(
                # Adding target releated fields:
                f.col("genes.*"),
                # Adding cell line specific annotation:
                f.array(f.col("diseaseCellLine")).alias("diseaseCellLines"),
                "biomarkerList",
                # Adding evidence specific stats:
                "phenotypicConsequenceLogFoldChange",
                "phenotypicConsequencePValue",
                "phenotypicConsequenceFDR",
                "statisticalMethod",
                "geneticInteractionPValue",
                "geneticInteractionScore",
                "geneInteractionType",
                # Adding disease information:
                f.lit(disease_from_source_mapped_id).alias("diseaseFromSourceMappedId"),
                f.lit(disease_from_source).alias("diseaseFromSource"),
                # Adding release releated fields:
                f.lit(self.release_version).alias("releaseVersion"),
                f.lit(self.release_date).alias("releaseDate"),
                # Data source specific fields:
                f.lit("ot_partner").alias("datatypeId"),
                f.lit("encore").alias("datasourceId"),
                f.lit("OTAR2062").alias("projectId"),
                f.lit("Encore project").alias("projectDescription"),
            )
            # Dropping all evidence for anchor and control genes:
            .filter(f.col("targetRole") == "library")
            .distinct()
        )

        # If there's a warning message for the sudy, add that:
        if "warningMessage" in parameters and parameters["warningMessage"] is not None:
            evidence_df = evidence_df.withColumn(
                "warningMessage", f.lit(parameters["warningMessage"])
            )

        return evidence_df


def main(outputFile: str, parameters: dict, cellPassportFile: str) -> None:

    # Initialize spark session
    spark = initialize_sparksession()

    # Opening and parsing the cell passport data from Sanger:
    cell_passport_df = GenerateDiseaseCellLines(cellPassportFile, spark).get_mapping()

    logging.info(f"Cell passport dataframe has {cell_passport_df.count()} rows.")
    logging.info("Parsing experiment data...")

    # Initialising evidence generator:
    evidenceGenerator = EncoreEvidenceGenerator(
        spark, cell_passport_df, parameters["sharedMetadata"]
    )

    # Create evidence for all experiments. Dataframes are collected in a list:
    evidence_dfs = []
    for experiment in parameters["experiments"]:
        evidence_dfs.append(evidenceGenerator.parse_experiment(experiment))

    # Filter out None values, so only dataframes with evidence are kept:
    evidence_dfs = list(filter(None, evidence_dfs))

    # combine all evidence dataframes into one:
    combined_evidence_df = reduce(lambda df1, df2: df1.union(df2), evidence_dfs)

    # The combined evidence is written to a file:
    write_evidence_strings(combined_evidence_df, outputFile)


if __name__ == "__main__":

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(
        description="This script parses ENCORE data and generates disease-target evidence."
    )
    parser.add_argument(
        "--output_file", "-o", type=str, help="Output file. gzipped JSON", required=True
    )
    parser.add_argument(
        "--parameter_file",
        "-p",
        type=str,
        help="A JSON file describing exeriment metadata",
        required=True,
    )
    parser.add_argument(
        "--data_folder",
        "-d",
        type=str,
        help="A folder in which the data files are located",
        required=True,
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="File into which the logs are saved",
        required=False,
    )
    parser.add_argument(
        "--cell_passport_file",
        "-c",
        type=str,
        help="File containing cell passport data",
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

    # Validating data folder:
    if not os.path.isdir(args.data_folder):
        logging.error(f"Data folder {args.data_folder} does not exist.")
        raise ValueError(f"Data folder {args.data_folder} does not exist.")

    # Reading parameter json:
    try:
        with open(args.parameter_file, "r") as parameter_file:
            parameters = json.load(parameter_file)
    except Exception as e:
        raise e(f"Could not read parameter file. {args.parameter_file}")

    # Updating parameters:
    parameters["sharedMetadata"]["data_folder"] = args.data_folder

    # Passing all the required arguments:
    main(args.output_file, parameters, args.cell_passport_file)
