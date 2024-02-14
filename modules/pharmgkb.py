#!/usr/bin/env python
"""This module adds a more granular description of the phenotype observed in the PharmGKB evidence."""

import argparse
import logging
import sys

import pyspark.sql.functions as f
from common.evidence import initialize_sparksession, write_evidence_strings
from common.ontology import add_efo_mapping
from pyspark.sql import SparkSession


def main(spark: SparkSession, pharmgkb_evidence_path: str, extracted_phenotypes_path: str, output_file_path: str, cache_dir: str) -> None:
    """This module overwrites the `phenotypeText` and adds `diseaseFromSourceMappedId` field in the PharmGKB evidence dataset.
    
    The original `phenotypeText` comes from PharmGKB directly, however it is tipically of little value for the user (more context in https://github.com/opentargets/curation/blob/master/docs/pharmacogenetics.md).

    Args:
        pharmgkb_evidence_path: Input gzipped JSON with the evidence submitted by ChEMBL.
        extracted_phenotypes_path: Input JSON with the phenotypes extracted from `genotypeAnnotationText`
        output_file_path: Output gzipped json file containing the pharmgkb evidence with a new `phenotypeText` and `diseaseFromSourceMappedId` fields.
        cache_dir: Directory to store the OnToma cache files in.
    """
    # Extract
    logging.info(f"PharmGKB evidence JSON file: {pharmgkb_evidence_path}")
    logging.info(f"Table of genotype descriptions to extracted phenotypes: {extracted_phenotypes_path}")
    spark.sparkContext.addFile(extracted_phenotypes_path)
    pharmgkb_df = spark.read.json(pharmgkb_evidence_path)
    pgx_phenotypes_df = spark.read.json(extracted_phenotypes_path.split("/")[-1])


    # Transform
    enriched_pharmgkb_df = (
        pharmgkb_df.drop("phenotypeText").join(pgx_phenotypes_df, on="genotypeAnnotationText", how="left")
        .withColumn("phenotypeText", f.explode("phenotypeText"))
        .persist()
    )
    enriched_pharmgkb_df = add_efo_mapping(
        enriched_pharmgkb_df.select("*", f.col("phenotypeText").alias("diseaseFromSource"), f.lit(None).alias("diseaseFromSourceId")),
        spark,
        cache_dir
    )
    logging.info("Disease mappings have been added.")

    # Load
    assert enriched_pharmgkb_df.count() >= pharmgkb_df.count(), "There are fewer evidence after processing."
    write_evidence_strings(enriched_pharmgkb_df, output_file_path)
    logging.info(f"{enriched_pharmgkb_df.count()} evidence strings have been saved to {output_file_path}. Exiting.")


def get_parser():
    """Get parser object for script pharmgkb.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--pharmgkb_evidence_path",
        help="Input gzipped JSON with the PharmGKB evidence submitted by EVA",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--extracted_phenotypes_path",
        help="Input TSV containing the categories of the clinical trial reason to stop. Instructions for applying the ML model here: https://github.com/ireneisdoomed/stopReasons.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_file_path", help="Output gzipped json file following the target safety liabilities data model.",
        type=str,
        required=True
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )
    return parser


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.StreamHandler(sys.stderr)

    spark = initialize_sparksession()
    args = get_parser().parse_args()

    main(
        spark,
        args.pharmgkb_evidence_path,
        args.extracted_phenotypes_path,
        args.output_file_path,
        args.cache_dir
    )
