#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the Broad CVDI Human Disease Portal."""

import argparse
import logging
import sys
from typing import Optional

import pandas as pd
import pyspark.sql.functions as f
from pandas import DataFrame as pd_dataframe
from pyspark.sql import DataFrame, SparkSession

from common.evidence import (
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping

METHOD_DESC = {
    "LOF + missense0.8 (MAF<0.1%)": "Mixed-effects test carried out with LOF and predicted-deleterious missense variants (missense score > 0.8) with a MAF smaller than 0.1%.",
    "LOF + missense0.5 (MAF<0.001%)": "Mixed-effects test carried out with LOF and predicted-deleterious missense variants (missense score > 0.5) with a MAF smaller than 0.001%.",
    "LOF (MAF<0.1%)": "Mixed-effects test carried out with LOF variants with a MAF smaller than 0.1%.",
    "LOF + missense0.5 (MAF<0.1%)": "Mixed-effects test carried out with LOF and predicted-deleterious missense variants (missense score > 0.5) with a MAF smaller than 0.1%.",
    "Cauchy": "Combined test after combining mask-specific using Cauchy distribution.",
}


def main(
    spark: SparkSession, input_path: str, ontoma_cache_dir: Optional[str] = None
) -> DataFrame:
    """
    This module extracts and processes target/disease evidence from the raw Broad CVDI Human Disease Portal.

    We use:
    - Table 6 as the main reference for significant target/disease associations. We filter out the ancestry specific associations because the study doesn't report any ancestry specific genes.
    - Table 3 contains the P cutoff for each of the methods. In the publication, they define statistical significance based on FDR < 0.01.

    Args:
        spark: Spark session
        input_path: Path to the input file
        ontoma_cache_dir: Path to the ontoma cache directory

    Returns:
        evd_df: DataFrame with CVDI's Human Disease Portal data following the t/d evidence schema
    """
    logging.info(f"File with the Broad CVDI Human Disease Portal data: {input_path}")

    # Load data
    associations_df = parse_associations_table(input_path)
    p_value_cutoff_df = parse_p_values_threshold_table(input_path)

    # Write output
    associations_df = (
        associations_df.merge(p_value_cutoff_df, left_on="method_name", right_on="Mask")
        .drop("Mask", axis=1)
        .drop_duplicates()
        # Dropping rows with no odds ratio:
        .astype({"OR [95%CI]": str})
    )
    evd_df = parse_evidence(spark, associations_df=associations_df)
    evd_df = add_efo_mapping(
        evidence_strings=evd_df, spark_instance=spark, ontoma_cache_dir=ontoma_cache_dir
    ).persist()

    if not 1_500 < evd_df.count() < 1_600:
        logging.exception(
            f"CVDI Disease Portal number of evidence are different from expected: {evd_df.count()}"
        )
        raise AssertionError(
            "CVDI Disease Portal number of evidence are different from expected."
        )
    logging.info(f"{evd_df.count()} evidence strings have been processed.")

    return evd_df


def parse_associations_table(input_path: str) -> pd_dataframe:
    """Parse Table 6 multiindex dataframe.

    Every column represents a different method

    We slice the dataframe to parse associations for each method and then merge them
    """

    def _slice_dataframe(df: pd_dataframe, method: str) -> pd_dataframe:
        """Slice a dataframe to extract the columns corresponding to a specific method.

        Args:
            df: DataFrame with columns indexed by method
            method: Method to extract

        Returns:
            df: DataFrame with columns corresponding to the specified method and a flatten structure of columns
        """
        df = df.xs(method, level=1, axis=1)
        df.columns = df.columns.get_level_values(1)
        return df.assign(method_name=method)

    associations_table = pd.read_excel(
        input_path,
        sheet_name="ST6",
        skiprows=1,
        header=[0, 1, 2],
        skipfooter=1,
    )[["phenotype", "Gene ID Ensembl", "Gene", "ALL ancestry"]]

    # Get the list of statistical models parsed in the column hierarchy
    statistical_models = list(
        {
            index_level_1
            for (index_level_1, _) in associations_table["ALL ancestry"].columns
        }
    )

    # Append index columns to each slice
    index_cols = ["phenotype", "Gene ID Ensembl", "Gene"]
    index_dataframe = associations_table[index_cols]
    index_dataframe.columns = index_dataframe.columns.get_level_values(0)
    return pd.concat(
        [
            pd.concat(
                [index_dataframe, _slice_dataframe(associations_table, model)], axis=1
            )
            for model in statistical_models
        ]
    )


def parse_p_values_threshold_table(input_path: str) -> pd_dataframe:
    """
    Parse Table 3 multiindex dataframe.

    Every column represents a different method
    """
    p_cutoff_table = pd.read_excel(
        input_path,
        sheet_name="ST3",
        skiprows=1,
        header=[0, 1],
        skipfooter=1,
    ).drop(["AoU", "UKB", "META (no correction)", "MGB"], axis=1)
    p_cutoff_table.columns = p_cutoff_table.columns.get_level_values(1)
    # Forward fill the 'Significance cutoff' values to fill the empty cells
    p_cutoff_table["Significance cutoff"] = p_cutoff_table[
        "Significance cutoff"
    ].ffill()
    return p_cutoff_table[p_cutoff_table["Significance cutoff"] == "FDR1%"].filter(
        ["Mask", "P cutoff"]
    )


def parse_evidence(spark: SparkSession, associations_df: pd_dataframe) -> DataFrame:
    """
    Parse CVDI's Human Disease Portal disease/target evidence.

    Args:
        associations_df: DataFrame with CVDI's Human Disease Portal associations already filtered and annotated with which p value cutoff is used per mask

    Returns:
        evd_df: DataFrame with CVDI's Human Disease Portal data following the t/d evidence schema (without EFO mapping)
    """
    return (
        spark.createDataFrame(
            associations_df,
        )
        .withColumn(
            "resourceScore",
            f.when(
                f.col("method_name") == f.lit("Cauchy"), f.col("Cauchy P-value")
            ).otherwise(f.col("Meta P-value")),
        )
        # Filter out non significant associations (thresholds vary depending on the mask)
        .filter(f.col("resourceScore") <= f.col("P cutoff"))
        .withColumn(
            "oddsRatio",
            f.split(f.regexp_extract(f.col("OR [95%CI]"), r"(\d+\.\d+)", 1), "\.")[
                0
            ].cast("float"),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            f.split(f.regexp_extract(f.col("OR [95%CI]"), r"\[(\d+\.\d+)", 1), "\.")[
                0
            ].cast("float"),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalUpper",
            f.split(f.regexp_extract(f.col("OR [95%CI]"), r"; (\d+\.\d+)\]", 1), "\.")[
                0
            ].cast("float"),
        )
        .withColumn("statisticalMethodOverview", f.col("method_name"))
        .replace(to_replace=METHOD_DESC, subset=["statisticalMethodOverview"])
        .select(
            f.lit("gene_burden").alias("datasourceId"),
            f.lit("genetic_association").alias("datatypeId"),
            f.lit("CVDI Human Disease Portal").alias("projectId"),
            f.lit("UK Biobank 450k/All of Us/MGB").alias("cohortId"),
            f.translate("phenotype", "_", " ").alias("diseaseFromSource"),
            f.lit(None).alias("diseaseFromSourceId"),
            f.col("Gene ID Ensembl").alias("targetFromSourceId"),
            f.col("cMAC").cast("int").alias("studyCasesWithQualifyingVariants"),
            f.lit(748879).alias("studySampleSize"),
            f.col("resourceScore"),
            f.col("oddsRatio"),
            f.col("oddsRatioConfidenceIntervalLower"),
            f.col("oddsRatioConfidenceIntervalUpper"),
            f.col("method_name").alias("statisticalMethod"),
            f.col("statisticalMethodOverview"),
            f.array(
                f.struct(
                    f.lit("Broad CVDI Human Disease Portal").alias("niceName"),
                    # E.g. https://hugeamp.org:8000/research.html?ancestry=mixed&cohort=UKB_450k_AoU_250k_MGB_53k_META_overlapcorrected&file=600Traits.csv&gene=MYBPC3&pageid=600_traits_app
                    f.concat(
                        f.lit(
                            "https://hugeamp.org:8000/research.html?ancestry=mixed&cohort=UKB_450k_AoU_250k_MGB_53k_META_overlapcorrected&file=600Traits.csv&gene="
                        ),
                        f.col("Gene"),
                        f.lit("&pageid=600_traits_app"),
                    ).alias("url"),
                )
            ).alias("urls"),
            f.array(f.lit("39210047")).alias("literature"),
        )
        .withColumn(
            "pValueExponent", f.log10(f.col("resourceScore")).cast("int") - f.lit(1)
        )
        .withColumn(
            "pValueMantissa",
            f.round(
                f.col("resourceScore") / f.pow(f.lit(10), f.col("pValueExponent")), 3
            ),
        )
        .withColumn(
            "studyCasesWithQualifyingVariants",
            f.when(
                f.col("studyCasesWithQualifyingVariants") == 0, f.lit(None)
            ).otherwise(f.col("studyCasesWithQualifyingVariants")),
        )
        .distinct()
    )


def get_parser():
    "Get parser object for script CvdiGeneBurden.py."
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input_path",
        help="Supplementary tables with CVDI's Human Disease Portal associations.",
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
        "--ontoma_cache_dir",
        help="Path to the ontoma cache directory.",
        type=str,
        required=False,
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
        input_path=args.input_path,
        ontoma_cache_dir=args.ontoma_cache_dir,
    )

    write_evidence_strings(evd_df, args.output)
    logging.info(f"Evidence strings have been saved to {args.output}. Exiting.")
