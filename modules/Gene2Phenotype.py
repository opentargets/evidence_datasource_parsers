#!/usr/bin/env python3
"""Evidence parser for Gene2Phenotype's disease panels."""

import argparse
import logging
import sys

from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession, functions as f, types as t

from common.ontology import add_efo_mapping
from common.evidence import initialize_sparksession, write_evidence_strings


G2P_mutationCsq2functionalCsq = {
    "uncertain": "SO_0002220",  # function_uncertain_variant
    "absent gene product": "SO_0002317",  # absent_gene_product
    "altered gene product structure": "SO_0002318",  # altered_gene_product_structure
    "5_prime or 3_prime UTR mutation": "SO_0001622",  # UTR_variant
    "increased gene product level": "SO_0002315",  # increased_gene_product_level
    "cis-regulatory or promotor mutation": "SO_0001566",  # regulatory_region_variant
}


def main(
    dd_file: str,
    eye_file: str,
    skin_file: str,
    cancer_file: str,
    cardiac_file: str,
    output_file: str,
    cache_dir: str,
    local: bool = False,
) -> None:
    # Initialize spark session
    spark = initialize_sparksession()

    # Read and process G2P's tables into evidence strings
    gene2phenotype_df = read_input_file(
        spark, dd_file, eye_file, skin_file, cancer_file, cardiac_file
    )
    logging.info(
        "Gene2Phenotype panels have been imported. Processing evidence strings."
    )

    evidence_df = process_gene2phenotype(gene2phenotype_df)

    evidence_df = add_efo_mapping(
        evidence_strings=evidence_df, ontoma_cache_dir=cache_dir, spark_instance=spark
    )
    logging.info("Disease mappings have been added.")

    # Saving data:
    write_evidence_strings(evidence_df, output_file)
    logging.info(
        f"{evidence_df.count()} evidence strings have been saved to {output_file}"
    )


def read_input_file(
    spark: SparkSession,
    dd_file: str,
    eye_file: str,
    skin_file: str,
    cancer_file: str,
    cardiac_file: str,
) -> DataFrame:
    """
    Reads G2P's panel CSV files into a Spark DataFrame forcing the schema
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
        .csv(
            [dd_file, eye_file, skin_file, cancer_file, cardiac_file],
            schema=gene2phenotype_schema,
            enforceSchema=True,
            header=True,
            sep=",",
            quote='"',
        )
    )


def process_gene2phenotype(gene2phenotype_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """

    return gene2phenotype_df.select(
        # Split pubmed IDs to list:
        f.split(f.col("pmids"), ";").alias("literature"),
        # Split phenotypes:
        f.split(f.col("phenotypes"), ";").alias("phenotypes"),
        # Renaming a few columns:
        f.col("gene symbol").alias("targetFromSourceId"),
        f.col("disease name").alias("diseaseFromSource"),
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
        f.regexp_replace(f.col("diseaseFromSource"), r".+-related ", "").alias(
            "diseaseFromSource"
        ),
        # Map functional consequences:
        translate(G2P_mutationCsq2functionalCsq)("mutation consequence").alias(
            "variantFunctionalConsequenceId"
        ),
        # Adding constant columns:
        f.lit("gene2phenotype").alias("datasourceId"),
        f.lit("genetic_literature").alias("datatypeId"),
    )


def translate(mapping):
    """
    Mapping consequences - to SO codes
    """

    def translate_(col):
        return mapping.get(col)

    return f.udf(translate_, t.StringType())


if __name__ == "__main__":
    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description="Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/"
    )
    parser.add_argument(
        "-d",
        "--dd_panel",
        help="DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-e",
        "--eye_panel",
        help="Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-s",
        "--skin_panel",
        help="Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-c",
        "--cancer_panel",
        help="Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-cr",
        "--cardiac_panel",
        help="Cardiac panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help="Absolute path of the gzipped, JSON evidence file.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-l", "--log_file", help="Filename to store the parser logs.", type=str
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )
    parser.add_argument(
        "--local",
        action="store_true",
        required=False,
        default=False,
        help="Flag to indicate if the script is executed locally or on the cluster",
    )

    args = parser.parse_args()

    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    cardiac_file = args.cardiac_panel
    output_file = args.output_file
    log_file = args.log_file
    cache_dir = args.cache_dir
    local = args.local

    # Configure logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if log_file:
        logging.config.fileConfig(fname=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Report input data:
    logging.info(f"DD panel file: {dd_file}")
    logging.info(f"Eye panel file: {eye_file}")
    logging.info(f"Skin panel file: {skin_file}")
    logging.info(f"Cancer panel file: {cancer_file}")
    logging.info(f"Cardiac panel file: {cardiac_file}")

    # Calling main:
    main(
        dd_file, eye_file, skin_file, cancer_file, cardiac_file, output_file, cache_dir
    )
