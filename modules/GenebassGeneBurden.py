#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the Genebass webservice."""

import argparse
import logging
import sys

import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame

from common.evidence import (
    import_trait_mappings,
    initialize_sparksession,
    write_evidence_strings,
)

METHOD_DESC = {
    "pLoF": "Burden test carried out with rare pLOF variants.",
    "missense|LC": "Burden test carried out with rare missense variants including low-confidence pLOF and in-frame indels.",
    "synonymous": "Burden test carried out with rare synonymous variants.",
    "pLoF|missense|LC": "Burden test carried out with pLOF or missense variants.",
}


def main(spark: SparkSession, genebass_data: str) -> DataFrame:
    """
    This module extracts and processes target/disease evidence from the raw Genebass Portal.
    """
    logging.info(f"File with the Genebass gene burden results: {genebass_data}")

    # Load data
    genebass_df = (
        SparkSession.getActiveSession()
        .read.parquet(genebass_data)
        .filter(F.col("Pvalue_Burden") <= 6.7e-7)
        .filter(F.col("trait_type") != "categorical")
        .select(
            "gene_id",
            "annotation",
            "n_cases",
            "n_controls",
            "trait_type",
            "phenocode",
            "description",
            "Pvalue_Burden",
            "BETA_Burden",
            "SE_Burden",
        )
        .distinct()
        .repartition(200)
        .persist()
    )

    # Write output
    evd_df = parse_genebass_evidence(spark, genebass_df)

    if evd_df.filter(F.col("resourceScore") == 0).count() != 0:
        logging.exception("There are evidence with a P value of 0.")
        raise AssertionError(
            f"There are {evd_df.filter(F.col('resourceScore') == 0).count()} evidence with a P value of 0."
        )
    if not 8_000 < evd_df.count() < 10_000:
        logging.exception(
            f"Genebass number of evidence are different from expected: {evd_df.count()} (expected between 8k and 10k)"
        )
        raise AssertionError("Genebass number of evidence are different from expected.")
    logging.info(f"{evd_df.count()} evidence strings have been processed.")

    return evd_df


def parse_genebass_evidence(spark: SparkSession, genebass_df: DataFrame) -> DataFrame:
    """
    Parse Genebass's disease/target evidence.
    Args:
        genebass_df: DataFrame with Genebass's portal data
    Returns:
        evd_df: DataFrame with Genebass's data following the t/d evidence schema.
    """
    to_keep = [
        "datasourceId",
        "datatypeId",
        "targetFromSourceId",
        "diseaseFromSource",
        "diseaseFromSourceId",
        "diseaseFromSourceMappedId",
        "pValueMantissa",
        "pValueExponent",
        "beta",
        "betaConfidenceIntervalLower",
        "betaConfidenceIntervalUpper",
        "resourceScore",
        "ancestry",
        "ancestryId",
        "projectId",
        "cohortId",
        "studySampleSize",
        "studyCases",
        "statisticalMethod",
        "statisticalMethodOverview",
        "literature",
    ]

    # WARNING: There are some associations with a p-value of 0.0 in Genebass.
    # This is a bug we still have to ellucidate and it might be due to a float overflow.
    # These evidence need to be manually corrected in order not to lose them and for them to pass validation
    # As an interim solution, their p value will equal to the minimum in the evidence set.
    logging.warning(
        f"There are {genebass_df.filter(F.col('Pvalue_Burden') == 0.0).count()} evidence with a p-value of 0.0."
    )
    minimum_pvalue = (
        genebass_df.filter(F.col("Pvalue_Burden") > 0.0)
        .agg({"Pvalue_Burden": "min"})
        .collect()[0]["min(Pvalue_Burden)"]
    )
    genebass_df = genebass_df.withColumn(
        "Pvalue_Burden",
        F.when(F.col("Pvalue_Burden") == 0.0, F.lit(minimum_pvalue)).otherwise(
            F.col("Pvalue_Burden")
        ),
    )

    return (
        genebass_df.withColumn("datasourceId", F.lit("gene_burden"))
        .withColumn("datatypeId", F.lit("genetic_association"))
        .withColumn("projectId", F.lit("Genebass"))
        .withColumn("cohortId", F.lit("UK Biobank 450k"))
        .withColumn("ancestry", F.lit("EUR"))
        .withColumn("ancestryId", F.lit("HANCESTRO_0009"))
        .withColumnRenamed("gene_id", "targetFromSourceId")
        .withColumnRenamed("description", "diseaseFromSource")
        .withColumnRenamed("phenocode", "diseaseFromSourceId")
        .join(
            import_trait_mappings(spark),
            on="diseaseFromSource",
            how="left",
        )
        .withColumnRenamed("Pvalue_Burden", "resourceScore")
        .withColumn(
            "pValueExponent",
            F.log10(F.col("resourceScore")).cast(T.IntegerType()) - F.lit(1),
        )
        .withColumn(
            "pValueMantissa",
            F.round(
                F.col("resourceScore") / F.pow(F.lit(10), F.col("pValueExponent")), 3
            ),
        )
        # Genebass reports effect sizes as beta, regardless of the type of measured trait, like other studies
        .withColumnRenamed("BETA_Burden", "beta")
        .withColumn("betaConfidenceIntervalLower", (F.col("beta") - F.col("SE_Burden")))
        .withColumn("betaConfidenceIntervalUpper", (F.col("beta") + F.col("SE_Burden")))
        .withColumn(
            "studySampleSize", (F.col("n_cases") + F.coalesce("n_controls", F.lit(0)))
        )
        .withColumnRenamed("n_cases", "studyCases")
        .withColumnRenamed("annotation", "statisticalMethod")
        .withColumn("statisticalMethodOverview", F.col("statisticalMethod"))
        .replace(to_replace=METHOD_DESC, subset=["statisticalMethodOverview"])
        .withColumn("literature", F.array(F.lit("36778668")))
        .select(to_keep)
        .distinct()
    )


def get_parser():
    "Get parser object for script GenebassGeneBurden.py."
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--genebass_data",
        help="Input parquet files with Genebass's burden associations.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Output gzipped json file following the gene_burden evidence data model.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--log_file",
        help="Destination of the logs generated by this script. Defaults to None",
        type=str,
        nargs="?",
        default=None,
    )

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    # Logger initializer. If no log_file is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    spark = initialize_sparksession()

    evd_df = main(
        spark=spark,
        genebass_data=args.genebass_data,
    )

    write_evidence_strings(evd_df, args.output)
    logging.info(f"Evidence strings have been saved to {args.output}. Exiting.")
