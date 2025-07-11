"""Evidence parser for ClinGen's Gene Validity Curations."""

import argparse
import logging

import pyspark.sql.functions as f
from pyspark.sql import DataFrame
from pyspark.sql.types import StringType, StructType, TimestampType

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping

logger = logging.getLogger(__name__)


def main(input_file: str, output_file: str, cache_dir: str) -> None:
    # Initialize spark session
    spark = initialize_sparksession()

    # Read and process Clingen's table into evidence strings

    clingen_df = read_input_file(input_file, spark_instance=spark)
    logger.info(
        "Gene Validity Curations table has been imported. Processing evidence strings."
    )

    evidence_df = process_clingen(clingen_df)

    evidence_df = add_efo_mapping(
        evidence_strings=evidence_df, spark_instance=spark, ontoma_cache_dir=cache_dir
    )
    logger.info("Disease mappings have been added.")

    write_evidence_strings(evidence_df, output_file)
    logger.info(
        f"{evidence_df.count()} evidence strings have been saved to {output_file}"
    )


def read_input_file(input_file: str, spark_instance) -> DataFrame:
    """
    Reads Gene Validity Curations CSV file into a Spark DataFrame forcing the schema
    The first 6 rows of this file include metadata that needs to be dropped
    """

    clingen_schema = (
        StructType()
        .add("GENE SYMBOL", StringType())
        .add("GENE ID (HGNC)", StringType())
        .add("DISEASE LABEL", StringType())
        .add("DISEASE ID (MONDO)", StringType())
        .add("MOI", StringType())
        .add("SOP", StringType())
        .add("CLASSIFICATION", StringType())
        .add("ONLINE REPORT", StringType())
        .add("CLASSIFICATION DATE", TimestampType())
        .add("GCEP", StringType())
    )

    clingen_df = (
        spark_instance.read.csv(input_file, schema=clingen_schema)
        # The first 6 rows of the GVC file include metadata that needs to be dropped
        .withColumn("idx", f.monotonically_increasing_id())
        .filter(f.col("idx") > 5)
        .drop("idx")
    )

    return clingen_df


def process_clingen(clingen_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """

    return (
        clingen_df.withColumn("datasourceId", f.lit("clingen"))
        .withColumn("datatypeId", f.lit("genetic_literature"))
        .withColumn("urls", f.struct(f.col("ONLINE REPORT").alias("url")))
        .select(
            "datasourceId",
            "datatypeId",
            f.trim(f.col("GENE SYMBOL")).alias("targetFromSourceId"),
            f.col("DISEASE LABEL").alias("diseaseFromSource"),
            f.col("DISEASE ID (MONDO)").alias("diseaseFromSourceId"),
            f.array(f.col("MOI")).alias("allelicRequirements"),
            f.array(f.col("urls")).alias("urls"),
            f.col("CLASSIFICATION").alias("confidence"),
            f.date_format(f.col("CLASSIFICATION DATE"), "yyyy-MM-dd").alias(
                "releaseDate"
            ),
            f.col("GCEP").alias("studyId"),
        )
    )


if __name__ == "__main__":
    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description="Parse ClinGen gene-disease associations from Gene Validity Curations"
    )
    parser.add_argument(
        "--input_file",
        help="Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_file",
        help="Absolute path of the gzipped JSON evidence file.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--log_file", type=str, help="Optional filename to redirect the logs into."
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )
    args = parser.parse_args()

    # Initialize logger:
    initialize_logger(__name__)

    # Report input data:
    logger.info(f"Clingen input file path: {args.input_file}")
    logger.info(f"Evidence output file path: {args.output_file}")
    logger.info(f"Cache directory: {args.cache_dir}")

    main(
        input_file=args.input_file,
        output_file=args.output_file,
        cache_dir=args.cache_dir,
    )
