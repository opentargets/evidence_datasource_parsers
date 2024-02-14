#!/usr/bin/env python3
"""Parser for fetch and format gene-essentiality dataset."""
from __future__ import annotations

import argparse
import logging
import logging.config
import sys
from typing import TYPE_CHECKING

from pyspark import SparkFiles
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t

from common.evidence import initialize_sparksession, write_evidence_strings

if TYPE_CHECKING:
    from pyspark.sql import Column


class DepMapEssentiality:
    """Parser for DepMap gene essentiality dataset."""

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
            .repartition("depmapId")
            # Joining with model information:
            .join(models_df, on="depmapId", how="left")
            # Joining expression data:
            .join(expression_df, on=["targetSymbol", "depmapId"], how="left")
            # Joining damaging mutation data:
            .join(damaging_mutation, on=["targetSymbol", "depmapId"], how="left")
            # Joining hotspot mutation data:
            .join(hotspots, on=["targetSymbol", "depmapId"], how="left")
            # Joining tissue mapping:
            .join(tissue_mapping, on="tissueFromSource", how="left")
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
                f.coalesce(f.col("tissueName"), f.lit("other")).alias("tissueName"),
                # Parsing mutation:
                f.when(f.col("damagingMutation").isNotNull(), "damaging")
                .when(f.col("hotspotMutation").isNotNull(), "hotspot")
                .alias("mutation"),
                # Measured values:
                "geneEffect",
                "expression",
                # Essentiality flag:
                f.coalesce(f.col("isEssential"), f.lit(False)).alias("isEssential"),
            )
            # Dropping rows with missing gene effect. This can happen when there's just no data for a gene in a given cell line:
            .filter(f.col("geneEffect").isNotNull())
            .persist()
        )

        # If only essential genes are needed, drop all other:
        if keep_only_essentials:
            self.essentials = self.essentials.filter(f.col("isEssential")).persist()

        # We have to assert the data is properly generated before returning:
        assert isinstance(self.essentials, DataFrame)

    def get_essential_genes_dataframe(self) -> DataFrame:
        """Returning aggregated view on essentiality.

        Returns:
            DataFrame: data grouped by gene then by tissue.

        Schema:
            root
            |-- targetSymbol: string (nullable = true)
            |-- isEssential: boolean (nullable = true)
            |-- depMapEssentiality: array (nullable = false)
            |    |-- element: struct (containsNull = false)
            |    |    |-- tissueId: string (nullable = true)
            |    |    |-- tissueName: string (nullable = true)
            |    |    |-- screens: array (nullable = false)
            |    |    |    |-- element: struct (containsNull = false)
            |    |    |    |    |-- depmapId: string (nullable = true)
            |    |    |    |    |-- cellLineName: string (nullable = true)
            |    |    |    |    |-- diseaseFromSource: string (nullable = true)
            |    |    |    |    |-- diseaseCellLineId: string (nullable = true)
            |    |    |    |    |-- mutation: string (nullable = true)
            |    |    |    |    |-- geneEffect: float (nullable = true)
            |    |    |    |    |-- expression: float (nullable = true)
        """
        # Return grouped essentiality data:
        return (
            # Aggregating data by gene:
            self.essentials.groupBy(
                "targetSymbol",
                "isEssential",
                "tissueId",
                "tissueName",
            )
            # Aggregating data further by tissue:
            .agg(
                f.collect_set(
                    f.struct(
                        f.col("depmapId").alias("depmapId"),
                        f.col("cellLineName").alias("cellLineName"),
                        f.col("diseaseFromSource").alias("diseaseFromSource"),
                        f.col("diseaseCellLineId").alias("diseaseCellLineId"),
                        f.col("mutation").alias("mutation"),
                        f.col("geneEffect").alias("geneEffect"),
                        f.col("expression").alias("expression"),
                    )
                ).alias("screens")
            )
            .groupBy("targetSymbol", "isEssential")
            .agg(
                f.collect_set(
                    f.struct(
                        f.col("tissueId").alias("tissueId"),
                        f.col("tissueName").alias("tissueName"),
                        f.col("screens").alias("screens"),
                    )
                ).alias("depMapEssentiality")
            )
            .persist()
        )

    def get_stats_on_essentials(self) -> None:
        """Print statistics on the essentiality dataset."""
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
            "file://" + SparkFiles.get(tissue_mapping_url.split("/")[-1]),
            sep=",",
            header=True,
        )

    def _read_and_melt(self, filename: str, value_name: str) -> DataFrame:
        """Read and melt the wide table into a long format.

        Damaging mutation and hotspot mutation tables are provided in a wide format, where each column represents a gene and each row is a cell line.
        The values indicate the number of mutations in the gene for the cell line.

        Args:
            filename (str): Path to the file.
            value_name (str): Name of the value column.

        Returns:
            DataFrame: Melted dataframe. Columns: depmapId, targetSymbol, value_name (e.g. geneEffect, expression, hotspotMutation, damagingMutation)
        """
        # Reading csv into dataframe:
        df = (
            self.spark.read.csv(filename, sep=",", header=True)
            # Some files don't have label for the first column:
            .withColumnRenamed("_c0", "depmapId")
            # Some files have the wrong label:
            .withColumnRenamed("ModelID", "depmapId")
        )

        # Extracting cell lines:
        genes = df.columns[1:]

        # Generate unpivot expression:
        unpivot_expression = f"""stack({len(genes)}, {', '.join([f"'{c}', `{c}`" for c in genes])}) as (gene_label, {value_name})"""

        # Transform dataset:
        return (
            df.select("depmapId", f.expr(unpivot_expression))
            .select(
                "depmapId",
                self._extract_gene_symbol(f.col("gene_label")),
                f.col(value_name).cast(t.FloatType()),
            )
            .repartition("depmapId")
        )

    @staticmethod
    def _extract_gene_symbol(gene_col: Column) -> Column:
        """Extract gene symbol from a string, where space separates gene symbol with other components.

        Data example:
            Essentials
            AAMP (14)
            AARS1 (16)

        Args:
            gene_col (Column): Column containing gene symbols.

        Returns:
            Column: gene symbol.
        """
        return f.split(gene_col, " ").getItem(0).alias("targetSymbol")

    def _prepare_models(self, model_file: str) -> DataFrame:
        """Prepare model data.

        Args:
            model_file (str): Path to the model file.

        Returns:
            DataFrame: columns: depmapId, cellLineName, modelId, tissueFromSource, diseaseFromSource
        """
        return self.spark.read.csv(model_file, sep=",", header=True).select(
            f.col("ModelID").alias("depmapId"),
            # If cell line name is provided, it's picked:
            f.when(f.col("CellLineName").isNotNull(), f.col("CellLineName"))
            # When not cell line name, but Cancer Cell Line Enciclopedia name is provided, that's picked:
            .when(f.col("CCLEName").isNotNull(), f.col("CCLEName"))
            # If none of these sources are available, the cell line name is generated from the disease name:
            .otherwise(f.concat(f.col("OncotreePrimaryDisease"), f.lit(" cells")))
            .alias("cellLineName"),
            f.col("SangerModelID").alias("modelId"),
            f.lower(f.col("OncotreeLineage")).alias("tissueFromSource"),
            f.col("OncotreePrimaryDisease").alias("diseaseFromSource"),
        )

    def _prepare_essential_genes(self, essential_gene_file: str) -> DataFrame:
        """Prepare essential gene list.

        This method reads the essential gene list table and returns a dataframe with gene symbols and essentiality flag:

        Daata example:
            Essentials
            AAMP (14)
            AARS1 (16)

        Args:
            essential_gene_file (str): Path to the essential gene list file.

        Returns:
            DataFrame: columns: targetSymbol, isEssential
        """
        return self.spark.read.csv(essential_gene_file, sep=",", header=True).select(
            self._extract_gene_symbol(f.col("Essentials")),
            f.lit(True).alias("isEssential"),
        )


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

    # Process essentiality from depmap:
    depmap_essentials = DepMapEssentiality(
        spark, depmap_input_folder, tissue_mapping_url, keep_essentials_only
    )

    # Write essentiality output:
    write_evidence_strings(
        depmap_essentials.get_essential_genes_dataframe(), output_file
    )


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
        "--log_file",
        help="Destination of the logs generated by this script. Defaults to None",
        type=str,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--essential_only",
        action="store_true",
        default=False,
        help="Only essential genes are kept in the output.",
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
    logging.info(f"Gene essentiality dataset is written to {args.output_file}.")
