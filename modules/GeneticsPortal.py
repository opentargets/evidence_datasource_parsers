#!/usr/bin/env python
"""This script pulls together data from Open Targets Genetics portal to generate disease/target evidence strings for the Platform."""

import argparse
import logging
import os
import sys
import tempfile

from psutil import virtual_memory
from pyspark.conf import SparkConf
from pyspark.sql import Column, DataFrame, SparkSession, Window
from pyspark.sql import functions as f
from pyspark.sql import types as t


def detect_spark_memory_limit():
    """Spark does not automatically use all available memory on a machine. When working on large datasets, this may
    cause Java heap space errors, even though there is plenty of RAM available. To fix this, we detect the total amount
    of physical memory and allow Spark to use (almost) all of it."""
    mem_gib = virtual_memory().total >> 30
    return int(mem_gib * 0.9)


def get_most_significant_qtl_effect(df: DataFrame) -> DataFrame:
    """Filter dataframe for rows with the most significant QTL effect.

    Args:
        df (DataFrame): Dataframe with the colocalization dataset

    Returns:
        DataFrame: From each group (defined by the groupbyColumns) only one row is returned.
    """
    # These are the columns we grouping by:
    group_columns = ["study_id", "gene_id", "chrom", "pos", "ref", "alt"]

    # Column name with the effect:
    qtl_significance_column = "qtlPValue"

    return (
        df
        # Ranking rows by the
        .withColumn(
            "effectRank",
            f.row_number().over(
                Window.partitionBy(*group_columns).orderBy(
                    f.col(qtl_significance_column).asc()
                )
            ),
        )
        # Fiter by rank:
        .filter(f.col("effectRank") == 1)
        # Drop helper columns:
        .drop("effectRank", "qtlPValue")
    )


def map_betas_to_so(beta_col: Column) -> Column:
    """Map effect values (betas) to sequence ontology (SO) terms.

    Args:
        beta_col (Column): occassionally missing float values

    Returns:
        Column: SO terms classification

    Examples:
        >>> columns = ['studyId', 'qtlEffect']
        >>> data = [('s1', -0.1),('s2', 12.0),('s3', None),('s4', 0.0)]
        >>> spark.createDataFrame(data,columns).withColumn('so', map_betas_to_so(f.col('qtlEffect'))).show()
        +-------+---------+----------+
        |studyId|qtlEffect|        so|
        +-------+---------+----------+
        |     s1|     -0.1|SO_0002316|
        |     s2|     12.0|SO_0002315|
        |     s3|     null|SO_0002314|
        |     s4|      0.0|SO_0002314|
        +-------+---------+----------+
        <BLANKLINE>
    """
    negative_beta = "SO_0002316"
    positive_beta = "SO_0002315"
    undefined_beta = "SO_0002314"

    return (
        f.when(beta_col < 0, negative_beta)
        .when(beta_col > 0, positive_beta)
        .otherwise(undefined_beta)
    )


def process_coloc(coloc_file: str) -> DataFrame:
    """For each GWAS loci, the direction of the biggest colocalizing QTL is returned.

    Args:
        coloc_file (str): parquet file pointing to the coloc dataset.

    Returns:
        DataFrame: gwas loci (study + variant id), qtl gene id + direction of effect
    """
    # "study_id", "gene_id", "chrom", "pos", "ref", "alt"
    # Filtering and processing coloc table:
    return (
        spark.read.parquet(coloc_file)
        .filter(
            # Dropping GWAS loci:
            (f.col("right_type") != "gwas")
            &
            # Excluding splice QTLs:
            (f.col("right_type") != "sqtl")
        )
        .select(
            f.col("left_chrom").alias("chrom"),
            f.col("left_pos").alias("pos"),
            f.col("left_ref").alias("ref"),
            f.col("left_alt").alias("alt"),
            f.col("left_study").alias("study_id"),
            f.col("right_gene_id").alias("gene_id"),
            f.col("left_var_right_study_beta").alias("qtlEffect"),
            f.col("left_var_right_study_pval").alias("qtlPValue"),
        )
        .distinct()
        # Windowing over the study/locus/gene QTLs and get the highest beta:
        .transform(get_most_significant_qtl_effect)
        .withColumn(
            "variantFunctionalConsequenceFromQtlId", map_betas_to_so(f.col("qtlEffect"))
        )
        .drop("qtlEffect")
        .persist()
    )


def get_parser():
    """Get parser object for script GeneticsPortal.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--locus2gene",
        help="Input table containing locus to gene scores.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--toploci",
        help="Table containing top loci for all studies.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--study", help="Table with all the studies.", type=str, required=True
    )
    parser.add_argument(
        "--variantIndex",
        help="Table with the variant indices (from gnomad 2.x).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--ecoCodes", help="Table with consequence ECO codes.", type=str, required=True
    )
    parser.add_argument(
        "--outputFile", help="Output gzipped json file.", type=str, required=True
    )
    parser.add_argument(
        "--threshold",
        help="Threshold applied on l2g score for filtering.",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--logFile",
        help="Destination of the logs generated by this script.",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--colocFile",
        help="Location of the coloc dataset in parquet format.",
        type=str,
        required=False,
    )

    return parser


def initialize_logger(logFile=None):
    """Logger initializer. If no logfile is specified, logs are written to stderr."""

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if logFile:
        logging.config.fileConfig(filename=logFile)
    else:
        logging.StreamHandler(sys.stderr)


def initialize_spark():
    """Spins up a Spark session."""

    # Initialize spark session
    spark_mem_limit = detect_spark_memory_limit()
    spark_conf = (
        SparkConf()
        .set("spark.driver.memory", f"{spark_mem_limit}g")
        .set("spark.executor.memory", f"{spark_mem_limit}g")
        .set("spark.driver.maxResultSize", "0")
        .set("spark.debug.maxToStringFields", "2000")
        .set("spark.sql.execution.arrow.maxRecordsPerBatch", "500000")
    )
    spark = (
        SparkSession.builder.config(conf=spark_conf).master("local[*]").getOrCreate()
    )
    logging.info(f"Spark version: {spark.version}")

    return spark


def load_eco_dict(vep_consequences: str):
    """
    Loads the csq to eco scores into a dict
    Returns: dict
    """

    # Load
    eco_df = spark.read.csv(
        vep_consequences, sep="\t", header=True, inferSchema=True
    ).select("Term", "Accession", f.col("eco_score").cast(t.DoubleType()))

    # Convert to python dict
    eco_dict = {}
    eco_link_dict = {}
    for row in eco_df.collect():
        eco_dict[row.Term] = row.eco_score
        eco_link_dict[row.Term] = row.Accession

    return (eco_dict, eco_link_dict)


def parse_genetics_evidence(genetics_df: DataFrame) -> DataFrame:
    """The JSON Schema format is applied to the df."""

    return (
        genetics_df.withColumn(
            "literature",
            f.when(
                f.col("pmid") != "",
                f.array(f.regexp_extract(f.col("pmid"), r"PMID:(\d+)$", 1)),
            )
            .when(f.col("study_id").contains("SAIGE"), f.array(f.lit("30104761")))
            .when(f.col("study_id").contains("FINNGEN"), f.array(f.lit("36653562"))),
        )
        .withColumn(
            "cohortId",
            f.when(
                f.col("study_id").contains("SAIGE"), f.array(f.lit("UK Biobank 500k"))
            ).when(
                f.col("study_id").contains("NEALE"), f.array(f.lit("UK Biobank 500k"))
            ),
        )
        .select(
            f.lit("ot_genetics_portal").alias("datasourceId"),
            f.lit("genetic_association").alias("datatypeId"),
            f.col("gene_id").alias("targetFromSourceId"),
            f.col("efo").alias("diseaseFromSourceMappedId"),
            f.col("literature"),
            f.col("pub_author").alias("publicationFirstAuthor"),
            "projectId",
            f.substring(f.col("pub_date"), 1, 4)
            .cast(t.IntegerType())
            .alias("publicationYear"),
            f.col("trait_reported").alias("diseaseFromSource"),
            f.col("study_id").alias("studyId"),
            f.col("sample_size").alias("studySampleSize"),
            f.col("pval_mantissa").alias("pValueMantissa"),
            f.col("pval_exponent").alias("pValueExponent"),
            f.col("odds_ratio").alias("oddsRatio"),
            f.col("oddsr_ci_lower").alias("oddsRatioConfidenceIntervalLower"),
            f.col("oddsr_ci_upper").alias("oddsRatioConfidenceIntervalUpper"),
            f.col("beta").alias("beta"),
            f.col("beta_ci_lower").alias("betaConfidenceIntervalLower"),
            f.col("beta_ci_upper").alias("betaConfidenceIntervalUpper"),
            f.col("y_proba_full_model").alias("resourceScore"),
            f.col("rsid").alias("variantRsId"),
            f.concat_ws(
                "_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")
            ).alias("variantId"),
            f.regexp_extract(f.col("consequence_link"), r"\/(SO.+)$", 1).alias(
                "variantFunctionalConsequenceId"
            ),
            "variantFunctionalConsequenceFromQtlId",
        )
        .dropDuplicates(
            ["variantId", "studyId", "targetFromSourceId", "diseaseFromSourceMappedId"]
        )
    )


def process_l2g_table(locus2gene: str, threshold: float) -> DataFrame:
    """Loads and processes locus-to-gene (L2G) score data."""

    return (
        spark.read.parquet(locus2gene)
        # Keep results trained on high or medium confidence gold-standards
        .filter(f.col("training_gs") == "high_medium")
        # Keep results from xgboost model
        .filter(f.col("training_clf") == "xgboost")
        # Keep rows with l2g score above the threshold:
        .filter(f.col("y_proba_full_model") >= threshold)
        # Only keep study, variant, gene and score info
        .select(
            "study_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gene_id",
            "y_proba_full_model",
        )
    )


def process_study_table(study_index: str) -> DataFrame:
    "Loads and processes disease information from the study table."

    return (
        spark.read.parquet(study_index)
        .select(
            "study_id",
            "pmid",
            "pub_date",
            "pub_author",
            "trait_reported",
            "trait_efos",
            f.col("n_initial").alias("sample_size"),
        )
        # Assign project based on the study author information
        .withColumn(
            "projectId",
            f.when(f.col("study_id").contains("FINNGEN"), "FINNGEN")
            .when(f.col("study_id").contains("NEALE"), "NEALE")
            .when(f.col("study_id").contains("SAIGE"), "SAIGE")
            .when(f.col("study_id").contains("GCST"), "GCST"),
        )
        # Warning! Not all studies have an EFO annotated (trait_efos is an empty array)
        # Also, some have multiple EFOs!
        # Studies with no EFO are kept, the array is exploded to capture each mapped trait
        .withColumn("efo", f.explode_outer(f.col("trait_efos")))
        .drop("trait_efos")
        # Drop records with HANCESTRO IDs as mapped trait
        .filter((~f.col("efo").contains("HANCESTRO")) | (f.col("efo").isNull()))
    )


def process_toploci_table(toploci: str) -> DataFrame:
    """Loads and processes association statistics (only pvalue is required) from top loci table."""

    return (
        spark.read.parquet(toploci)
        .select(
            "study_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "beta",
            "beta_ci_lower",
            "beta_ci_upper",
            "pval_mantissa",
            "pval_exponent",
            "odds_ratio",
            "oddsr_ci_lower",
            "oddsr_ci_upper",
        )
        # Problem: Large OR values which cannot be represented with a single precision float are
        # automatically casted to 'infinity' when the df is exported to JSON, hence failing validation.
        # This was also a problem for ES, as these evidence were not being loaded (https://github.com/opentargets/platform/issues/1687)
        # Decision: OR is set to null.
        .withColumn(
            "odds_ratio",
            f.when(
                f.col("odds_ratio").between(sys.float_info.min, sys.float_info.max),
                f.col("odds_ratio"),
            ),
        )
        .withColumn(
            "oddsr_ci_lower",
            f.when(
                f.col("oddsr_ci_lower").between(sys.float_info.min, sys.float_info.max),
                f.col("oddsr_ci_lower"),
            ),
        )
        .withColumn(
            "oddsr_ci_upper",
            f.when(
                f.col("oddsr_ci_upper").between(sys.float_info.min, sys.float_info.max),
                f.col("oddsr_ci_upper"),
            ),
        )
    )


def process_consequences_table(variant_index: str, vep_consequences: str) -> DataFrame:
    """Loads and processes variant's most severe functional consequence (SO ID)."""

    eco_dicts = spark.sparkContext.broadcast(load_eco_dict(vep_consequences))

    get_consequence_link_udf = f.udf(lambda x: eco_dicts.value[1][x], t.StringType())

    get_most_severe_consequence_udf = f.udf(
        # Extract most sereve csq per gene.
        # Create UDF that reverse sorts csq terms using eco score dict, then select
        # the first item. Then apply UDF to all rows in the data.
        lambda arr: sorted(
            arr, key=lambda x: eco_dicts.value[0].get(x, 0), reverse=True
        )[0],
        t.StringType(),
    )

    return (
        spark.read.parquet(variant_index)
        # Explode consequences, only keeping canonical transcript
        .selectExpr(
            "chrom_b38 as chrom",
            "pos_b38 as pos",
            "ref",
            "alt",
            "vep.most_severe_consequence as most_severe_csq",
            """explode(
                filter(vep.transcript_consequences, x -> x.canonical == 1)
            ) as tc
            """,
        )
        # Keep required fields from consequences struct
        .selectExpr(
            "chrom",
            "pos",
            "ref",
            "alt",
            "most_severe_csq",
            "tc.gene_id as gene_id",
            "tc.consequence_terms as csq_arr",
        )
        # Get most severe consequences
        .withColumn(
            "most_severe_gene_csq", get_most_severe_consequence_udf(f.col("csq_arr"))
        )
        .withColumn(
            "consequence_link", get_consequence_link_udf(f.col("most_severe_gene_csq"))
        )
    )


def process_variant_rsid(variant_index: str):
    """Load and extract rsIDs and genomic coordinates from the variant index."""

    return (
        spark.read.parquet(variant_index)
        # chrom_b38|pos_b38
        # Explode consequences, only keeping canonical transcript
        .selectExpr("chrom_b38 as chrom", "pos_b38 as pos", "ref", "alt", "rsid")
    )


def write_evidence_strings(evidence_df: DataFrame, output_file: str) -> None:
    """Exports the table to a compressed JSON file containing the evidence strings."""
    (
        evidence_df.toPandas().to_json(
            output_file,
            orient="records",
            lines=True,
            compression="gzip",
        )
    )


def main(
    locus2gene: str,
    toploci: str,
    study_index: str,
    variant_index: str,
    vep_consequences: str,
    threshold: float,
    output_file: str,
    coloc_table: str,
):
    logging.info(f"Locus2gene table: {locus2gene}")
    logging.info(f"Coloc table: {coloc_table}")
    logging.info(f"Top locus table: {toploci}")
    logging.info(f"Study table: {study_index}")
    logging.info(f"Variant index table: {variant_index}")
    logging.info(f"ECO code table: {vep_consequences}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"l2g score threshold: {threshold}")

    # Load and process the input files into dataframes
    l2g_df = process_l2g_table(locus2gene, threshold)
    pvals_df = process_toploci_table(toploci)
    studies_df = process_study_table(study_index)
    variant_consequences_df = process_consequences_table(
        variant_index, vep_consequences
    )
    variant_rsid_df = process_variant_rsid(variant_index)
    coloc_df = process_coloc(coloc_table)

    # Join datasets together
    genetics_df = (
        l2g_df
        # Join L2G to pvals, using study and variant info as key
        .join(pvals_df, on=["study_id", "chrom", "pos", "ref", "alt"])
        # Join this to the study info, using study_id as key
        .join(studies_df, on="study_id", how="inner")
        # Join transcript consequences
        .join(
            variant_consequences_df,
            on=["chrom", "pos", "ref", "alt", "gene_id"],
            how="left",
        )
        # Bring rsIDs
        .join(variant_rsid_df, on=["chrom", "pos", "ref", "alt"], how="left")
        # Filling missing consequences
        .fillna(
            {
                "most_severe_gene_csq": "intergenic_variant",
                "consequence_link": "http://purl.obolibrary.org/obo/SO_0001628",
            }
        )
        # Joining with colocalizing QTL effects:
        .join(
            coloc_df,
            on=["study_id", "gene_id", "chrom", "pos", "ref", "alt"],
            how="left",
        )
    )

    # Write output
    logging.info("Evidence strings have been processed. Saving...")
    genetics_df = parse_genetics_evidence(genetics_df)
    write_evidence_strings(genetics_df, output_file)
    return 0


if __name__ == "__main__":
    args = get_parser().parse_args()
    initialize_logger(args.logFile)

    global spark
    spark = initialize_spark()

    main(
        locus2gene=args.locus2gene,
        toploci=args.toploci,
        study_index=args.study,
        variant_index=args.variantIndex,
        vep_consequences=args.ecoCodes,
        threshold=args.threshold,
        output_file=args.outputFile,
        coloc_table=args.colocFile,
    )
