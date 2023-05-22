import argparse
import logging
from typing import List

import pandas as pd
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import DataFrame

from common.evidence import initialize_sparksession, write_evidence_strings

PROBES_SETS = [
    "Bromodomains chemical toolbox",
    "Chemical Probes for Understudied Kinases",
    "Chemical Probes.org",
    "Gray Laboratory Probes",
    "High-quality chemical probes",
    "MLP Probes",
    "Nature Chemical Biology Probes",
    "Open Science Probes",
    "opnMe Portal",
    "Probe Miner (suitable probes)",
    "Protein methyltransferases chemical toolbox",
    "SGC Probes",
    "Tool Compound Set",
    "Concise Guide to Pharmacology 2019/20",
    "Kinase Chemogenomic Set (KCGS)",
    "Kinase Inhibitors (best-in-class)",
    "Novartis Chemogenetic Library (NIBR MoA Box)",
    "Nuisance compounds in cellular assays",
]


def collapse_cols_data_in_array(
    df: DataFrame, source_cols: List[str], destination_col: str
) -> DataFrame:
    """Collapses the data in a single column when the information is one-hot encoded.

    Args:
        df (DataFrame): Dataframe containing the data for the different probes
        source_cols (List[str]): List of columns containing the data to be collapsed
        destination_col (str): Name of the column where the array will be stored

    Returns:
        DataFrame: Dataframe with a new column containing the sources that have data for a specific probe

    Examples:
        >>> probes_df = spark.createDataFrame([("A", 1, 0), ("B", 1, 1)], ["probe", "source1", "source2"])
        >>> collapse_cols_data_in_array(df, ["source1", "source2"], "datasource").show()
        +-----+-------+-------+------------------+
        |probe|source1|source2|        datasource|
        +-----+-------+-------+------------------+
        |    A|      1|      0|         [source1]|
        |    B|      1|      1|[source1, source2]|
        +-----+-------+-------+------------------+
        <BLANKLINE>

    """
    # Escape the name of the columns in case they contain spaces
    source_cols = [f"`{e}`" for e in source_cols if " " in e]
    return df.withColumn(
        destination_col,
        f.array([f.when(df[c] == 1, c.replace(r"`", "")) for c in source_cols]),
    ).withColumn(
        destination_col, f.array_except(f.col(destination_col), f.array(f.lit(None)))
    )


def clean_origin_col():
    """Removes the substring ' probe' from the origin column to just state if the probe has been reported from an experimental or computational approach."""
    return f.array_distinct(
        f.expr("transform(origin, x -> trim(regexp_replace(x, ' probe', '')))")
    )


def extract_hq_flag():
    """Returns a flag indicating if the probe is high-quality or not."""
    return f.when(
        f.array_contains(f.col("datasourceIds"), "High-quality chemical probes"),
        f.lit(True),
    ).otherwise(f.lit(False))


def convert_stringified_array_to_array(col_name: str) -> DataFrame:
    """Converts a column of stringified arrays to an array column.

    Args:
        col_name: Name of the column that contains the stringified array

    Examples:
        >>> df = spark.createDataFrame([("['BI-1829', 'No control']")], t.StringType())
        >>> df.select(convert_stringified_array_to_array("value").alias("res")).show(truncate=False)
        +---------------------+
        |res                  |
        +---------------------+
        |[BI-1829, No control]|
        +---------------------+
        <BLANKLINE>

    """
    return f.split(f.translate(col_name, "[]'", ""), ", ").cast(
        t.ArrayType(t.StringType())
    )


def replace_dash(col_name):
    """Converts to null those values that only contain `-`."""
    return f.when(f.col(col_name).cast(t.StringType()) != "-", f.col(col_name))


def process_probes_data(probes_excel: str) -> List[DataFrame]:
    """Metadata about the compound and the scores given by the different sources."""
    return (
        spark.createDataFrame(
            pd.read_excel(probes_excel, sheet_name="PROBES", header=0, index_col=0)
            # Probes that do not have an associated target are marked as nulls
            .query("target.notnull()").reset_index()
        )
        # Collect list of datasources for each probe
        .transform(
            lambda df: collapse_cols_data_in_array(df, PROBES_SETS, "datasourceIds")
        )
        # Collecting the list of detection methods of the probe
        .transform(
            lambda df: collapse_cols_data_in_array(
                df,
                ["experimental probe", "calculated probe"],
                "origin",
            )
        ).select(
            "pdid",
            f.col("compound_name").alias("id"),
            clean_origin_col().alias("origin"),
            # Flag the high-quality probes and remove this from the list of datasources
            extract_hq_flag().alias("isHighQuality"),
            f.explode(
                f.array_except(
                    f.col("datasourceIds"),
                    f.array(f.lit("High-quality chemical probes")),
                )
            ).alias("datasourceId"),
        )
    )


def process_probes_targets_data(probes_excel: str) -> DataFrame:
    """Collection of targets associated with the probes and their scores."""
    return (
        spark.createDataFrame(
            pd.read_excel(
                probes_excel, sheet_name="PROBES TARGETS", header=0, index_col=0
            )
            # Probes that do not have an associated target are marked with "-"
            .query("gene_name != '-'").reset_index()
        )
        .filter(f.col("organism") == "Homo sapiens")
        .withColumn(
            "mechanismOfAction",
            f.when(
                f.col("action") != "-",
                f.split(f.col("action"), ";"),
            ),
        )
        .select(
            "pdid",
            f.col("target").alias("targetFromSource"),
            f.col("inchikey").alias("inchiKey"),
            "mechanismOfAction",
            replace_dash("`P&D probe-likeness score`")
            .alias("probesDrugsScore")
            .cast(t.DoubleType()),
            replace_dash("`Probe Miner Score`")
            .alias("probeMinerScore")
            .cast(t.DoubleType()),
            replace_dash("`Cells score (Chemical Probes.org)`")
            .alias("scoreInCells")
            .cast(t.DoubleType()),
            replace_dash("`Organisms score (Chemical Probes.org)`")
            .alias("scoreInOrganisms")
            .cast(t.DoubleType()),
        )
    )


def process_drug_index(drug_index: str) -> DataFrame:
    """OT Drug index to resolve the compound names to a ChEMBLID."""
    return spark.read.parquet(drug_index).selectExpr("id as drugId", "inchiKey")


def process_probes_sets_data(probes_excel: str) -> DataFrame:
    """Metadata about the different sources of probes."""
    return (
        spark.createDataFrame(
            pd.read_excel(
                probes_excel, sheet_name="COMPOUNDSETS", header=0, index_col=0
            )
        )
        .selectExpr("COMPOUNDSET as datasourceId", "SOURCE_URL as url")
        .filter(f.col("url").startswith("http"))
    )


def process_targets_data(probes_excel: str) -> DataFrame:
    """Look-up table between the gene symbols and the UniProt IDs."""
    return spark.createDataFrame(
        pd.read_excel(
            probes_excel, sheet_name="TARGETS", header=0, index_col=0
        ).reset_index()
    ).selectExpr("target as targetFromSource", "uniprot as targetFromSourceId")


def process_controls_data(controls_lut: str) -> DataFrame:
    """Look-up table between the probes and their control compounds.

    Note: We have asked P&Ds to include this info in the main probes
    Excel so that we don't have to process their SQL dump.
    """

    return (
        spark.read.csv(controls_lut, header=True)
        .select(
            f.col("probeid").alias("pdid"),
            convert_stringified_array_to_array("control").alias("control"),
        )
        # Many cases have "1" in the control array
        # This indicates the probe's efficacy has been tested with a control
        # But we don't know which one
        .withColumn(
            "control",
            f.when(~f.array_contains("control", f.lit("1")), f.col("control")),
        )
        .withColumn("control", f.explode_outer("control"))
    )


def main(probes_excel: str, drug_index: str, probes_controls: str) -> DataFrame:
    """Main logic of the script."""

    probes_df = process_probes_data(probes_excel)
    probes_targets_df = process_probes_targets_data(probes_excel)
    probes_sets_df = process_probes_sets_data(probes_excel)
    probes_controls_df = process_controls_data(probes_controls)
    targets_df = process_targets_data(probes_excel)
    drug_idx = process_drug_index(drug_index)

    grouping_cols = [
        "targetFromSourceId",
        "id",
        "drugId",
        "mechanismOfAction",
        "origin",
        "control",
        "isHighQuality",
        "probesDrugsScore",
        "probeMinerScore",
        "scoreInCells",
        "scoreInOrganisms",
    ]

    return (
        probes_targets_df.join(probes_df, on="pdid", how="left")
        .join(drug_idx, on="inchiKey", how="left")
        .join(targets_df, on="targetFromSource", how="left")
        .join(probes_sets_df, on="datasourceId", how="left")
        .join(probes_controls_df, on="pdid", how="left")
        .groupBy(grouping_cols)
        .agg(
            f.collect_set(
                f.struct(
                    f.col("datasourceId").alias("niceName"), f.col("url").alias("url")
                )
            ).alias("urls")
        )
    )


def get_parser():
    """Get parser object for script ChemicalProbes.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--probes_excel_path",
        help="Path to the Probes&Drugs Probes XLSX that contains the main tables.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--drug_index",
        help="Directory of parquet files that stores OT's disease index.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--probes_controls",
        help="CSV that contains the export of the `probe` SQL table to extract info about the probes' controls.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Output gzipped json file following the chemical probes data model.",
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

    logging.basicConfig(
        level=logging.INFO,
        format="%(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    global spark
    spark = initialize_sparksession()

    out = main(
        probes_excel=args.probes_excel_path,
        probes_controls=args.probes_controls,
        drug_index=args.drug_index,
    )
    write_evidence_strings(out, args.output)
    logging.info(f"Probes dataset has been saved to {args.output}. Exiting.")
