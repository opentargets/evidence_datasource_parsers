"""Ingests GTEx V8 data and computes the baseline expression data."""

import argparse
import functools
import gzip
import json

import pandas as pd

from . import metrics


def read_gtex_data(gtex_source_data, low_expression_threshold):
    """Read and preprocess GTEx source data."""
    # Read the source data. It's in the GCT format, which is a regular TSV, but the first two lines are metadata which we don't need.
    df = pd.read_csv(gtex_source_data, sep="\t", skiprows=2)
    # Drop column with gene names, which we don't need at all.
    df.drop(columns=["Description"], inplace=True)
    # Drop Ensembl gene version: "ENSG00000223972.2" → "ENSG00000223972".
    df["Name"] = df["Name"].str.split(".").str[0]
    # Ensembl gene IDs are used as index throughout the entire processing code.
    df.set_index(keys=["Name"], inplace=True)
    return df


def filter_out_low_expression_genes(df, threshold):
    """Genes with very low expression tend to generate noisy data. One particular case is an expression vector where one tissue has expression of 0.1 TPM and all the other ones are 0.0 exactly. This causes the gene to be called as extremely specifically expressed, even though most of these issues are just false discoveries.
    A common threshold used for filtering out genes with very low expression is to filter out those where expression in all tissues is less than 1.0 TPM. This is also used by HPA, and supported by internal investigations into quality control.
    """
    max_expr_across_tissues = df.max(axis=1)
    expr_level_mask = max_expr_across_tissues >= threshold
    return df[expr_level_mask]


def calculate_specificity_metrics(df, threshold):
    a = df.iloc[:, :0].copy()
    a["gini"] = df.apply(metrics.gini_coefficient, axis=1)
    a["hpaSpecificity"] = df.apply(
        functools.partial(metrics.hpa_specificity, low_expression_threshold=threshold), axis=1
    )
    a["hpaDistribution"] = df.apply(
        functools.partial(metrics.hpa_distribution, low_expression_threshold=threshold), axis=1
    )


def pack_adatiss_row(row):
    """Given a row with Adatiss results for a given gene, pack them into the list of (tissue, z-score) values ready for output."""
    # Extract the column names and values from the row
    cols = adatiss.columns
    vals = row.values.tolist()

    # Pack the values into a list of dictionaries
    dicts = []
    for col, val in zip(cols, vals):
        dicts.append(
            {
                "bodyPartLevel": "tissue",
                "bodyPartName": col,
                "bodyPartId": name_mapping[col],
                "adatissScore": round(val, 3),
            }
        )

    return dicts


def remove_adatiss_if_none(d):
    if d["adatissScores"] != d["adatissScores"]:  # using the fact that NaN != NaN
        return {k: v for k, v in d.items() if k != "adatissScores"}
        print("removed")
    else:
        return d


def main(gtex_source_data, low_expression_threshold):
    df = read_gtex_data(gtex_source_data)
    df = filter_out_low_expression_genes(df, low_expression_threshold)
    a = calculate_specificity_metrics(df, low_expression_threshold)

    # Calculate expression specificity scores.

    # Pack the data for output.
    df_name_mapping = pd.read_table(
        "/home/ktsukanov/repositories/open-targets/curation/mappings/biosystem/gtex_v8_tissues.tsv"
    )
    name_mapping = {
        gtex_tissue: ontology_code
        for gtex_tissue, ontology_code in df_name_mapping[["gtex_tissue_name", "ontology_code"]].values.tolist()
    }
    # Verify that all columns are present in the mapping file.
    for tissue in df.columns:
        if tissue not in name_mapping:
            raise AssertionError(f"GTEx tissue {tissue} is not present in the provided mappings file")

    # Map annotation names.
    annotation_column_map = {
        "gini": "gini",
        "hpa_specificity": "hpaSpecificity",
        "hpa_distribution": "hpaDistribution",
    }
    columns_to_filter = list(annotation_column_map.keys())
    a_out = a[columns_to_filter].rename(columns=annotation_column_map)
    a_out["gini"] = a_out["gini"].round(3)
    a_out["adatissScores"] = adatiss.apply(pack_adatiss_row, axis=1)

    # Create a list to store JSON objects
    json_list = []

    # Iterate over each row index
    for idx in df.index:
        # Create a dictionary for the current row
        row_dict = {
            "ensemblGeneId": idx,
            "expression": [
                {
                    "bodyPartLevel": "tissue",
                    "bodyPartId": key,
                    "bodyPartName": name_mapping[key],
                    "tpm": value,
                }
                for key, value in df.loc[idx].to_dict().items()
            ],
            "expressionSpecificity": remove_adatiss_if_none(a_out.loc[idx].to_dict()),
        }
        # Append the dictionary to the list
        json_list.append(row_dict)

    # Output the data.
    with gzip.open("expression_v2.1.jsonl.gz", "wt", compresslevel=9) as f:
        for obj in json_list:
            f.write(json.dumps(obj))
            f.write("\n")


parser = argparse.ArgumentParser()
parser.add_argument("--gtex-source-data", required=True, type=str, help="GTEx by-gene median TPM counts, in GCT format")
parser.add_argument(
    "--low-expression-threshold",
    required=True,
    type=float,
    help="Genes which are not expressed at least this much in at least one tissue will be filtered out",
)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.gtex_source_data, args.low_expression_threshold)
