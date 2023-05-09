#!/usr/bin/env python3
"""Parser for Target Enabling Packages (TEP) downloaded from the Structural Genomics Consortium website."""
from __future__ import annotations
from typing import TYPE_CHECKING


import argparse
import logging
import logging.config
import sys
from functools import reduce
from pyspark.sql import SparkSession, functions as f, types as t, DataFrame, Column

from common.evidence import write_evidence_strings, initialize_sparksession
from pyspark.sql import types as t, functions as f
from pyspark import SparkFiles

if TYPE_CHECKING:
    from pyspark.sql import DataFrame, SparkSession


class DepMapEssentiality:

    # List of files to be used:
    MODELS = "Model.csv"
    ESSENTIAL_GENE_LIST = "CRISPRInferredCommonEssentials.csv"
    GENE_EFFECTS = "CRISPRGeneEffect.csv"
    GENE_EXPRESSION = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    MUTATION_HOTSPOTS = "OmicsSomaticMutationsMatrixHotspot.csv"
    MUTATION_DAMAGING = "OmicsSomaticMutationsMatrixDamaging.csv"

    def __init__(
        self: DepMapEssentiality,
        spark: SparkSession,
        input_folder: str,
        tissue_mapping_url: str,
        keep_only_essentials: bool = True,
    ) -> None:

        self.spark = spark

        # Reading model data:
        models_df = self._prepare_models(f"{input_folder}/{self.MODELS}")

        # Get a list of essential genes:
        essential_genes_df = self._prepare_essential_genes(
            f"{input_folder}/{self.ESSENTIAL_GENE_LIST}"
        )

        # Get a list of effect for each gene/cell line pair:
        gene_effect_df = self._read_and_melt(
            f"{input_folder}/{self.GENE_EFFECTS}", "geneEffect"
        )

        # Read and process expression data:
        expression_df = self._read_and_melt(
            f"{input_folder}/{self.GENE_EXPRESSION}", "expression"
        )

        # Read and process hotspot mutation:
        hotspots = self._read_and_melt(
            f"{input_folder}/{self.MUTATION_HOTSPOTS}", "hotspotMutation"
        ).filter(f.col("hotspotMutation") > 0)

        # Read and process damaging mutation:
        damaging_mutation = self._read_and_melt(
            f"{input_folder}/{self.MUTATION_DAMAGING}", "damagingMutation"
        ).filter(f.col("damagingMutation") > 0)

        # Get tissue mapping:
        tissue_mapping = self._get_tissue_mapping(tissue_mapping_url)

        # Join all the data together:
        self.essentials = (
            # The melted effect table:
            gene_effect_df
            # Narrowing down the gene effect table to genes that are in the essential gene list:
            .join(essential_genes_df, on="targetSymbol", how="left")
            # Joining with model information:
            # .join(models_df, on="depmapId", how="left")
            .join(models_df, on="depmapId", how="inner")
            # Joining expression data:
            .join(expression_df, on=["targetSymbol", "depmapId"], how="left")
            # Joining damaging mutation data:
            .join(damaging_mutation, on=["targetSymbol", "depmapId"], how="left")
            # Joining hotspot mutation data:
            .join(hotspots, on=["targetSymbol", "depmapId"], how="left")
            # Joining tissue mapping:
            .join(tissue_mapping, on="oncotreeLineage", how="left")
            # Format data:
            .select(
                "targetSymbol",
                "depmapId",
                # Cell-line fields:
                "cellLineName",
                f.col("modelId").alias("diseaseCellLineId"),
                # Disease:
                "diseaseFromSource",
                # Tissue fields:
                "tissueId",
                f.when(f.col("tissueName").isNull(), "other")
                .otherwise(f.col("tissueName"))
                .alias("tissueName"),
                # Parsing mutation:
                f.when(f.col("damagingMutation").isNotNull(), "damaging")
                .when(f.col("hotspotMutation").isNotNull(), "hotspot")
                .alias("mutation"),
                # Measured values:
                "geneEffect",
                "expression",
                # Essentiality flag:
                f.when(f.col("isEssential") == True, True)
                .otherwise(False)
                .alias("isEssential"),
            ).persist()
        )

        # If only essential genes are needed, drop all other:
        if keep_only_essentials:
            self.essentials = self.essentials.filter(
                f.col("isEssential") == True
            ).persist()

        # We have to assert the data is properly generated before returning:
        assert isinstance(self.essentials, DataFrame)

    def get_essential_genes_dataframe(self) -> DataFrame:
        return self.essentials

    def get_stats_on_essentials(self) -> None:
        entry_count = self.essentials.count()
        gene_count = self.essentials.select("targetSymbol").distinct().count()
        disease_count = self.essentials.select("diseaseFromSource").distinct().count()

        print(f"Number of entries: {entry_count}")
        print(f"Number of essential genes: {gene_count}")
        print(f"Number of unique diseases: {disease_count}")

    def _get_tissue_mapping(self, tissue_mapping_url: str) -> DataFrame:
        """Fetch tissue mapping stored as a tsv in github URL.

        Args:
            tissue_mapping_url (str): URL to the tissue mapping file.

        Returns:
            DataFrame: columns:
        """
        self.spark.sparkContext.addFile(tissue_mapping_url)
        return self.spark.read.csv(
            SparkFiles.get(tissue_mapping_url.split("/")[-1]), sep=",", header=True
        )

    def _read_and_melt(self, filename: str, value_name: str) -> DataFrame:
        # Reading csv into dataframe:
        df = self.spark.read.csv(filename, sep=",", header=True).withColumnRenamed(
            "_c0", "depmapId"
        )

        # Extracting cell lines:
        genes = df.columns[1:]

        # Generate unpivot expression:
        unpivot_expression = f"""stack({len(genes)}, {', '.join([f"'{c}', `{c}`" for c in genes])}) as (gene_label, {value_name})"""

        # Transform dataset:
        return df.select("depmapId", f.expr(unpivot_expression)).select(
            "depmapId",
            self._extract_gene_symbol(f.col("gene_label")),
            f.col(value_name).cast(t.FloatType()),
        )

    @staticmethod
    def _extract_gene_symbol(gene_col: Column) -> Column:
        return f.split(gene_col, " ").getItem(0).alias("targetSymbol")

    def _prepare_models(self, model_file: str) -> DataFrame:
        return self.spark.read.csv(model_file, sep=",", header=True).select(
            f.col("ModelID").alias("depmapId"),
            f.col("CellLineName").alias("cellLineName"),
            f.col("SangerModelID").alias("modelId"),
            f.col("OncotreeLineage").alias("oncotreeLineage"),
            f.col("OncotreePrimaryDisease").alias("diseaseFromSource"),
        )

    def _prepare_essential_genes(self, essential_gene_file: str) -> DataFrame:
        return self.spark.read.csv(essential_gene_file, sep=",", header=True).select(
            self._extract_gene_symbol(f.col("Essentials")),
            f.lit(True).alias("isEssential"),
        )


def get_depmap_essentials(
    spark: SparkSession,
    input_folder: str,
    tissue_mapping_url: str,
    keep_essentials_only: bool = False,
) -> DataFrame:
    """Wrapper function around the DepMap gene essentiality parser.

    Args:
        spark (SparkSession):
        input_folder (str): input folder
        tissue_mapping_url (str): path to tissue mapping
        keep_essentials_only (bool): flag indicating if only essential genes are needed.

    Returns:
        DataFrame: Parsed essentiality data in spark dataframe.
    """
    depmap_essentials = DepMapEssentiality(
        spark, input_folder, tissue_mapping_url, keep_essentials_only
    )
    depmap_essentials.get_stats_on_essentials()
    return depmap_essentials.get_essential_genes_dataframe()


def main(
    depmap_input_folder: str,
    tissue_mapping_url: str,
    output_file: str,
    keep_essentials_only: bool = False,
) -> None:
    """Processing gene essentiality evidene based on different sources.

    Args:
        depmap_input_folder (str): Path to the input files for DepMap
        tissue_mapping_url (str): Path to the tissue mappings URL
        output_file (str): path for the output file (gzipped json)
        keep_essentials_only (bool): flag indicating if only essential genes are needed.
    """

    logging.info(f"DepMap input folder: {depmap_input_folder}")
    logging.info(f"Tissue mapping URL: {tissue_mapping_url}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"Keep only essential genes: {keep_essentials_only}")

    # Initializing spark session:
    spark = initialize_sparksession()

    gene_essentiality = [
        # Generate gene essentiality tables from depmap data:
        get_depmap_essentials(
            spark, depmap_input_folder, tissue_mapping_url, keep_essentials_only
        ).persist(),
        # Potential further sources will come here:
    ]

    # merging all dataset into a single table:
    gene_essentiality_df = reduce(
        lambda df1, df2: df1.unionByName(df2, allowMissingColumns=True),
        gene_essentiality,
    )

    # Write output:
    write_evidence_strings(gene_essentiality_df, output_file)


def parse_command_line_parameters():
    parser = argparse.ArgumentParser(
        description="Generating gene essentiality annotation."
    )

    parser.add_argument(
        "--depmap_input_folder",
        help="Path to folder in which input files are stored from DepMap.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--depmap_tissue_mapping",
        help="Cell/tissue to UBERON mappings for DepMap.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_file",
        help="Gene essentiality annotation in compressed json.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--essential_only",
        help="Only essential genes are kept in the output.",
        type=bool,
        required=False,
        default=False,
    )
    parser.add_argument(
        "--log_file",
        help="Destination of the logs generated by this script. Defaults to None",
        type=str,
        default=None,
        required=False,
    )

    return parser


if __name__ == "__main__":

    # Reading output file name from the command line:
    args = parse_command_line_parameters().parse_args()

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if args.log_file:
        logging.Handler().setFormatter(
            fmt="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s"
        )
        logging.config.fileConfig(args.log_file, disable_existing_loggers=False)
    else:
        logging.StreamHandler(sys.stderr)

    # Generating gene essentiality dataset:
    main(
        args.depmap_input_folder,
        args.depmap_tissue_mapping,
        args.output_file,
        args.essential_only,
    )
