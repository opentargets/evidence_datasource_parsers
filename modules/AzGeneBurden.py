#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the raw data published in PMID:34662886."""

import argparse
import logging
import sys

from pyspark.sql.dataframe import DataFrame
from pyspark.sql.functions import array, format_string, col, lit, split, when
from pyspark.sql.types import DoubleType, IntegerType

from common.evidence import initialize_sparksession, write_evidence_strings

METHOD_DESC = {
    "ptvraredmg": "Burden test carried out with singleton pLOF variants.",
    "UR": "Burden test carried out with pLOFs with a MAF smaller than 0.001%.",
    "raredmg": "Burden test carried out with pLOFs with a MAF smaller than 0.01%.",
    "URmtr": "Burden test carried out with pLOFs with a MAF smaller than 0.1%.",
    "ptv": "Burden test carried out with pLOFs with a MAF smaller than 1%.",
    "rec": "Burden test carried out with singleton pLOFs and deleterious missense variants.",
    "syn": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.001%.",
    "raredmgmtr": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.01%.",
    "ptv5pcnt": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.1%.",
    "flexdmg": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 1%.",
    "flexnonsynmtr": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 1%.",
}


def main(az_binary_data: str, az_quant_data: str, output_file: str) -> None:
    """
    This module extracts and processes target/disease evidence from the raw data published in PMID:34662886.
    Args:
    """
    logging.info(f"File with the AZ PheWAS Portal binary traits associations: {az_binary_data}")
    logging.info(f"File with the AZ PheWAS Portal quantitative traits associations: {az_quant_data}")

    # Load data
    az_phewas_df = (
        spark.read.parquet(az_binary_data)
        .withColumnRenamed("BinOddsRatioLCI", "LCI")
        .withColumnRenamed("BinOddsRatioUCI", "UCI")
        .unionByName(spark.read.parquet(az_quant_data), allowMissingColumns=True)
        .filter(col('pValue') <= 2e-9)
        .distinct()
        .repartition(20)
        .persist()
    )

    assert 12000 < az_phewas_df.count() < 13000, "AZ PheWAS Portal number of evidence are different from expected."

    # Write output
    evd_df = parse_az_phewas_evidence(az_phewas_df)
    assert evd_df.filter(col('resourceScore') == 0).count() >= 0, "P-value is 0 for some associations."
    logging.info('Evidence strings have been processed. Saving...')

    write_evidence_strings(evd_df, output_file)

    logging.info(f"{evd_df.count()} evidence strings have been saved to {output_file}. Exiting.")


def parse_az_phewas_evidence(az_phewas_df: DataFrame) -> DataFrame:
    """
    Parse Astra Zeneca's PheWAS Portal evidence.
    Args:
        az_phewas_df: DataFrame with the Regeneron data
    Returns:
        evd_df: DataFrame with the Regeneron data following the t/d evidence schema.
    """
    TO_KEEP = [
        "datasourceId",
        "datatypeId",
        "targetFromSourceId",
        "diseaseFromSource",
        "diseaseFromSourceMappedId",
        "pValueMantissa",
        "pValueExponent",
        "beta",
        "betaConfidenceIntervalLower",
        "betaConfidenceIntervalUpper",
        "oddsRatio",
        "oddsRatioConfidenceIntervalLower",
        "oddsRatioConfidenceIntervalUpper",
        "resourceScore",
        "ancestry",
        "ancestryId",
        "literature",
        "projectId",
        "cohortId",
        "studySampleSize",
        "studyCases",
        "statisticalMethod",
        "statisticalMethodOverview",
        "urls"
    ]

    BASE_URL = "https://azphewas.com/geneView/7e2a7fab-97f0-45f7-9297-f976f7e667c8/%s/glr/%s"

    return (
        az_phewas_df.withColumn("datasourceId", lit("gene_burden"))
        .withColumn("datatypeId", lit("genetic_association"))
        .withColumn("literature", array(lit("34375979")))
        .withColumn("projectId", lit("AstraZeneca PheWAS Portal"))
        .withColumn("cohortId", lit("UK Biobank"))
        .withColumnRenamed("Gene", "targetFromSourceId")
        .withColumnRenamed("Phenotype", "diseaseFromSource")
        .withColumn("diseaseFromSourceMappedId", lit("placeholder"))
        .withColumn("resourceScore", col("pValue").cast(DoubleType()))
        .withColumn("pValueMantissa", split(col("pValue"), "e").getItem(0).cast(DoubleType()))
        .withColumn("pValueExponent", split(col("pValue"), "e").getItem(1).cast(IntegerType()))
        .withColumn(
            "beta",
            when(col("Type") == "Quantitative", col("beta")),
        )
        .withColumn(
            "betaConfidenceIntervalLower",
            when(col("Type") == "Quantitative", col("LCI")),
        )
        .withColumn(
            "betaConfidenceIntervalUpper",
            when(col("Type") == "Quantitative", col("UCI")),
        )
        .withColumn(
            "oddsRatio",
            when(col("Type") == "Binary", col("binOddsRatio")),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            when(col("Type") == "Binary", col("LCI")),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalUpper",
            when(col("Type") == "Binary", col("UCI")),
        )
        .withColumn("ancestry", lit("EUR"))
        .withColumn("ancestryId", lit("HANCESTRO_0005"))
        .withColumnRenamed("nSamples", "studySampleSize")
        .withColumn(
            "studyCases", when(col("Type") == "Quantitative", col("studySampleSize")).otherwise(col("BinNcases"))
        )
        .withColumnRenamed("CollapsingModel", "statisticalMethod")
        .withColumn("statisticalMethodOverview", col("statisticalMethod"))
        .replace(to_replace=METHOD_DESC, subset=['statisticalMethodOverview'])
        .withColumn(
            "urls",
            when(col("Type") == "Binary", format_string(BASE_URL, col("targetFromSourceId"), lit("binary"))).otherwise(
                format_string(BASE_URL, col("targetFromSourceId"), lit("continuous"))
            ),
        )
        .select(TO_KEEP)
        .distinct()
    )


def get_parser():
    """Get parser object for script RegeneronGeneBurden.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--az_binary_data',
        help='Input parquet files with AZ\'s PheWAS associations of binary traits.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--az_quant_data',
        help='Input parquet files with AZ\'s PheWAS associations of quantitative traits.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Output gzipped json file following the target safety liabilities data model.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--log_file',
        help='Destination of the logs generated by this script. Defaults to None',
        type=str,
        nargs='?',
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

    global spark
    spark = initialize_sparksession()

    main(az_binary_data=args.az_binary_data, az_quant_data=args.az_quant_data, output_file=args.output)
