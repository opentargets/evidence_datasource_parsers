#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the raw data published in PMID:34662886."""

import argparse
import logging
import sys

from numpy import nan
from pandas import read_excel
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame
from pyspark.sql.functions import aggregate, array, element_at, expr, col, concat, lit, split, translate, udf, when
from pyspark.sql.types import ArrayType, DoubleType, IntegerType

from common.evidence import get_exponent, get_mantissa, initialize_sparksession, write_evidence_strings

ANCESTRY_TO_ID = {"EUR": "HANCESTRO_0005", "EAS": "HANCESTRO_0009", "AFR": "HANCESTRO_0010", "SAS": "HANCESTRO_0006"}

MARKER_TO_METHOD_DESC = {
    "M1.singleton": "Burden test carried out with singleton pLOF variants.",
    "M1.0001": "Burden test carried out with pLOFs with a MAF smaller than 0.001%.",
    "M1.001": "Burden test carried out with pLOFs with a MAF smaller than 0.01%.",
    "M1.01": "Burden test carried out with pLOFs with a MAF smaller than 0.1%.",
    "M1.1": "Burden test carried out with pLOFs with a MAF smaller than 1%.",
    "M3.singleton": "Burden test carried out with singleton pLOFs and deleterious missense variants.",
    "M3.0001": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.001%.",
    "M3.001": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.01%.",
    "M3.01": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 0.1%.",
    "M3.1": "Burden test carried out with pLOFs and deleterious missense with a MAF smaller than 1%.",
}


def main(regeneron_data: str, gwas_studies: str, spark_instance: SparkSession) -> DataFrame:
    """
    This module extracts and processes target/disease evidence from the raw data published in PMID:34662886.
    Args:
    """
    logging.info(f"File with the results from the Regeneron burden analyses: {regeneron_data}")
    logging.info(f"File with GWAS Catalog studies: {gwas_studies}")

    # Load data

    TO_KEEP = [
        "Gene",
        "Trait",
        "Trait type",
        "UKB detailed trait name",
        "Trait_String",
        "Marker",
        "Marker type",
        "P-value",
        "Effect (95% CI)",
        "N cases with 0|1|2 copies of effect allele",
        "N controls with 0|1|2 copies of effect allele",
        "Ancestry",
    ]

    eur_evd_df = (
        spark_instance.createDataFrame(
            read_excel(regeneron_data, sheet_name="SD2").filter(items=TO_KEEP).replace({nan: None}).astype(str)
        )
        .withColumn("Ancestry", lit("EUR"))
        .withColumnRenamed("UKB detailed trait name", "Study tag")
    )

    non_eur_evd_df = spark_instance.createDataFrame(
        read_excel(regeneron_data, sheet_name="SD3").filter(items=TO_KEEP).replace({nan: None}).astype(str)
    ).withColumnRenamed("Trait_String", "Study tag")

    summary_stats_df = spark_instance.createDataFrame(
        read_excel(regeneron_data, sheet_name="SD4", skiprows=0, header=1)
        .iloc[2:]
        .reset_index(drop=True)
        .filter(items=["Study tag", "Study Accession"])
        .replace({nan: None})
        .astype(str)
    ).withColumn("Study tag", split("Study tag", "__SV").getItem(0))

    gwas_studies_df = (
        spark_instance.read.csv(gwas_studies, header=True, inferSchema=True, sep="\t")
        .filter(col("LINK").contains("34662886"))
        .select(
            col("STUDY ACCESSION").alias("Study Accession"),
            element_at(split("MAPPED_TRAIT_URI", "/"), -1).alias("MAPPED_TRAIT"),
        )
    )
    assert gwas_studies_df.count() > 7900, "Downloaded GWAS studies data are not complete."

    regeneron_df = (
        # Combine European with Non European associations
        eur_evd_df.union(non_eur_evd_df)
        # Keep only gene based analyses
        .filter(expr("`Marker type` == 'Burden'"))
        # add study accession from sumstats
        .join(summary_stats_df, on="Study tag", how="left")
        # add mapped trait from GWAS Catalog
        .join(gwas_studies_df, on="Study Accession", how="left")
        .filter(col('P-Value') <= 2.18e-11)
        .distinct()
        .persist()
    )

    evd_df = parse_regeneron_evidence(regeneron_df)
    if evd_df.filter(col('resourceScore') == 0).count() != 0:
        logging.error("There are evidence with a P value of 0.")
        raise AssertionError(
            f"There are {evd_df.filter(col('resourceScore') == 0).count()} evidence with a P value of 0."
        )

    if not 7000 < evd_df.count() < 8000:
        logging.error(f"REGENERON number of evidence are different from expected: {evd_df.count()}")
        raise AssertionError("AZ PheWAS Portal number of evidence are different from expected.")

    if evd_df.filter(col('studySampleSize').isNull()).count() != 0:
        logging.error(f"There are {evd_df.filter(col('studySampleSize').isNull()).count()} evidence with a null sample size.")
        raise AssertionError(
            f"Study sample size is missing for some associations."
        )
    logging.info(f"{evd_df.count()} evidence strings have been processed.")

    return evd_df


def parse_regeneron_evidence(regeneron_df: DataFrame) -> DataFrame:
    """
    Parse REGENERON's ExWAS analyses.
    Args:
        regeneron_df: DataFrame with the Regeneron data
    Returns:
        evd_df: DataFrame with the Regeneron data following the t/d evidence schema.
    """
    to_keep = [
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
        "studyId",
        "studySampleSize",
        "studyCases",
        "statisticalMethod",
        "statisticalMethodOverview",
    ]

    get_exponent_udf = udf(get_exponent, IntegerType())
    get_mantissa_udf = udf(get_mantissa, DoubleType())

    return (
        regeneron_df.withColumn("datasourceId", lit("gene_burden"))
        .withColumn("datatypeId", lit("genetic_association"))
        .withColumn("literature", array(lit("34662886")))
        .withColumn("projectId", lit("REGENERON"))
        .withColumn("cohortId", lit("UK Biobank"))
        .withColumnRenamed("Gene", "targetFromSourceId")
        .withColumnRenamed("Trait", "diseaseFromSource")
        .withColumnRenamed("MAPPED_TRAIT", "diseaseFromSourceMappedId")
        .withColumn("resourceScore", col("P-value").cast(DoubleType()))
        .withColumn("pValueMantissa", get_mantissa_udf(col("P-value")))
        .withColumn("pValueExponent", get_exponent_udf(col("P-value")))
        # Parse interval by removing unwanted characters and splitting.
        # Ex: "-0.097 (-0.125, -0.069)" -> [-0.097, -0.125, -0.069]
        .withColumn("effectParsed", split(translate(col('Effect (95% CI)'), '(),', ''), ' '))
        .withColumn(
            "beta",
            when(
                col("Trait type") == "QT",
                col("effectParsed").getItem(0).cast(DoubleType()),
            ),
        )
        .withColumn(
            "betaConfidenceIntervalLower",
            when(
                col("Trait type") == "QT",
                col("effectParsed").getItem(1).cast(DoubleType()),
            ),
        )
        .withColumn(
            "betaConfidenceIntervalUpper",
            when(
                col("Trait type") == "QT",
                col("effectParsed").getItem(2).cast(DoubleType()),
            ),
        )
        .withColumn(
            "oddsRatio",
            when(
                col("Trait type") == "BT",
                col("effectParsed").getItem(0).cast(DoubleType()),
            ),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            when(
                col("Trait type") == "BT",
                col("effectParsed").getItem(1).cast(DoubleType()),
            ),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalUpper",
            when(
                col("Trait type") == "BT",
                col("effectParsed").getItem(2).cast(DoubleType()),
            ),
        )
        .withColumnRenamed("ancestry", "Ancestry")
        .withColumn("ancestryId", col("ancestry"))
        .replace(to_replace=ANCESTRY_TO_ID, subset=['ancestryId'])
        .withColumnRenamed("Study Accession", "studyId")
        # Parse number of individuals and sum them up. Bear in mind that studyControls can be null
        # Ex: "421865|1618|2" -> 423485
        # Ex: "NA|NA|NA" -> null
        .withColumn(
            "studyCases",
            aggregate(
                split(col('N cases with 0|1|2 copies of effect allele'), '\|').cast(ArrayType(IntegerType())),
                lit(0),
                lambda acc, x: acc + x,
            ),
        )
        .withColumn(
            "studyControls",
            aggregate(
                split(col('N controls with 0|1|2 copies of effect allele'), '\|').cast(ArrayType(IntegerType())),
                lit(0),
                lambda acc, x: acc + x,
            ),
        )
        .withColumn(
            "studySampleSize",
            when(col('studyControls').isNotNull(), col('studyCases') + col('studyControls')).otherwise(
                col('studyCases')
            ),
        )
        .withColumn("statisticalMethod", concat(lit("ADD-WGR-FIRTH_"), col("Marker")))
        .withColumn("statisticalMethodOverview", col("Marker"))
        .replace(to_replace=MARKER_TO_METHOD_DESC, subset=['statisticalMethodOverview'])
        .select(to_keep)
        .distinct()
    )


def get_parser():
    """Get parser object for script RegeneronGeneBurden.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--regeneron_data',
        help='Input Excel file with the data published in PMID:34662886.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--gwas_studies',
        help='Input TSV containing all GWAS Catalog studies.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Output gzipped json file following the gene_burden evidence data model.',
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

    spark = initialize_sparksession()

    evd_df = main(regeneron_data=args.regeneron_data, gwas_studies=args.gwas_studies, spark_instance=spark)

    write_evidence_strings(evd_df, args.output)
    logging.info(f"Evidence strings have been saved to {args.output}. Exiting.")
