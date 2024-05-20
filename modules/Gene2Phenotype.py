#!/usr/bin/env python3
"""Evidence parser for Gene2Phenotype's disease panels."""

import argparse
import logging
from collections import OrderedDict
from typing import List, Optional

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping

G2P_mutationCsq2functionalCsq = OrderedDict(
    [
        ("uncertain", "SO_0002220"),
        ("cis-regulatory or promotor mutation", "SO_0001566"),
        ("5_prime or 3_prime UTR mutation", "SO_0001622"),
        ("increased gene product level", "SO_0002315"),
        ("decreased gene product level", "SO_0002316"),
        ("altered gene product structure", "SO_0002318"),
        ("absent gene product", "SO_0002317"),
    ]
)


def main(
    gene2phenotype_panels: List[str],
    output_file: str,
    cache_dir: Optional[str],
) -> None:
    # Get logger:
    logger = logging.getLogger(__name__)

    # Initialize spark session
    spark = initialize_sparksession()

    # Log the processed panels:
    for panel in gene2phenotype_panels:
        logger.info(f"Processing the following Gene2Phenotype panel: {panel}")

    # Read and process G2P's tables into evidence strings
    gene2phenotype_df = read_input_files(spark, gene2phenotype_panels)
    logger.info("Processing evidence strings.")

    evidence_df = process_gene2phenotype(gene2phenotype_df)

    evidence_df = add_efo_mapping(
        evidence_strings=evidence_df, ontoma_cache_dir=cache_dir, spark_instance=spark
    )
    logger.info("Disease mappings have been added.")

    # Saving data:
    write_evidence_strings(evidence_df, output_file)
    logger.info(
        f"{evidence_df.count()} evidence strings have been saved to {output_file}"
    )


def read_input_files(
    spark: SparkSession,
    gene2phenotype_panels: List[str],
) -> DataFrame:
    """Read G2P's panel CSV files into a Spark DataFrame forcing the schema.

    Args:
        spark: The Spark session.
        gene2phenotype_panels (List[str]): List of paths to G2P's panel CSV files.

    Returns:
        DataFrame: The DataFrame with the parsed data.
    """

    gene2phenotype_schema = (
        t.StructType()
        .add("gene symbol", t.StringType())
        .add("gene mim", t.IntegerType())
        .add("disease name", t.StringType())
        .add("disease mim", t.StringType())
        .add("confidence category", t.StringType())
        .add("allelic requirement", t.StringType())
        .add("mutation consequence", t.StringType())
        .add("phenotypes", t.StringType())
        .add("organ specificity list", t.StringType())
        .add("pmids", t.StringType())
        .add("panel", t.StringType())
        .add("prev symbols", t.StringType())
        .add("hgnc id", t.IntegerType())
        .add("gene disease pair entry date", t.TimestampType())
        .add("cross cutting modifier", t.StringType())
        .add("mutation consequence flag", t.StringType())
        .add("confidence value flag", t.StringType())
        .add("comments", t.StringType())
        .add("variant consequence", t.StringType())
        .add("disease ontology", t.StringType())
    )

    return (
        spark.read.option("multiLine", True)
        .option("encoding", "UTF-8")
        # There are rows with newlines and quotes, so there must be some escaping:
        .option("quote", '"')
        .option("escape", '"')
        .csv(
            gene2phenotype_panels,
            schema=gene2phenotype_schema,
            enforceSchema=True,
            header=True,
            sep=",",
            multiLine=True,
        )
    )


def process_gene2phenotype(gene2phenotype_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """
    return gene2phenotype_df.select(
        # Split pubmed IDs to list:
        f.split(f.col("pmids"), ";").alias("literature"),
        # Renaming a few columns:
        f.col("gene symbol").alias("targetFromSourceId"),
        f.col("panel").alias("studyId"),
        f.col("confidence category").alias("confidence"),
        # Parsing allelic requirements:
        f.when(
            f.col("allelic requirement").isNotNull(),
            f.array(f.col("allelic requirement")),
        ).alias("allelicRequirements"),
        # Parsing disease from source identifier:
        f.when(f.col("disease ontology").isNotNull(), f.col("disease ontology"))
        .when(
            (~f.col("disease mim").contains("No disease mim"))
            & (f.col("disease mim").isNotNull()),
            f.concat(f.lit("OMIM:"), f.col("disease mim")),
        )
        .otherwise(f.lit(None))
        .alias("diseaseFromSourceId"),
        # Cleaning disease names:
        f.regexp_replace(f.col("disease name"), r".+-related ", "").alias(
            "diseaseFromSource"
        ),
        # Map functional consequences:
        parse_functional_consequence(G2P_mutationCsq2functionalCsq)(
            "mutation consequence"
        ).alias("variantFunctionalConsequenceId"),
        # Adding constant columns:
        f.lit("gene2phenotype").alias("datasourceId"),
        f.lit("genetic_literature").alias("datatypeId"),
    )


def parse_functional_consequence(mapping):
    """Map semicolon separated list of consequence terms to SO codes.

    Args:
        mapping(Dict[str]): ordered dictionary with consequence mappings.

    Return:
        udf
    """

    def __translate(col):
        return sorted(
            # Get a list of consequence labels mapped to SO terms:
            [
                mapping.get(consequence_term)
                for consequence_term in col.split(";")
                if consequence_term in mapping
            ],
            # Sort SO terms based on the order in the consequnce map:
            key=lambda x: list(G2P_mutationCsq2functionalCsq.values()).index(x),
            # Get the most severe SO term:
        )[-1]

    return f.udf(__translate, t.StringType())


def parse_parameters() -> argparse.Namespace:
    """Parse CLI arguments."""
    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description="Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/"
    )
    parser.add_argument(
        "--panels",
        help="Space separated list of gene2phenotype panels downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "--output_file",
        help="Absolute path of the gzipped, JSON evidence file.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--log_file", help="Filename to store the parser logs.", type=str
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    # Get parameters:
    args = parse_parameters()

    gene2phenotype_panels: List[str] = args.panels
    output_file: str = args.output_file
    log_file: Optional[str] = args.log_file
    cache_dir: Optional[str] = args.cache_dir

    # Configure logger:
    initialize_logger(__name__, log_file)

    # Calling main:
    main(
        gene2phenotype_panels,
        output_file,
        cache_dir,
    )
