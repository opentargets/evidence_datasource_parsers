"""Collection of functions to calculate various expression specificity metrics."""

import functools
import os
import subprocess

import numpy as np
import pandas as pd


# General expression metrics.


def gini_coefficient(row):
    """Calculate the Gini coefficient of a Pandas row. Adapted from
    https://github.com/oliviaguest/gini/blob/master/gini.py"""
    assert min(row) >= 0
    if max(row) == 0:
        return np.nan
    # Values must be sorted.
    array = np.sort(row)
    # Index per array element.
    index = np.arange(1, array.shape[0] + 1)
    # Number of array elements.
    n = array.shape[0]
    # Gini coefficient.
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))


def hpa_specificity(row, low_expression_threshold):
    """HPA specificity metric. See: https://www.proteinatlas.org/about/assays+annotation#classification_rna."""
    expr = sorted(row)
    if expr[-1] < low_expression_threshold:
        return "Not detected"
    if expr[-2] == 0 or expr[-1] / expr[-2] >= 4.0:
        return "Tissue enriched"
    for i in range(2, 6):
        if expr[-i - 1] == 0 or expr[-i] / expr[-i - 1] >= 4.0:
            return "Group enriched"
    mean = sum(expr) / len(row)
    if 1 <= sum([e / mean >= 4.0 for e in expr]) <= 5:
        return "Tissue enhanced"
    return "Low tissue specificity"


def hpa_distribution(row, low_expression_threshold):
    """HPA distribution metric. See: https://www.proteinatlas.org/about/assays+annotation#classification_rna."""
    expr = sorted(row)
    if expr[-1] < low_expression_threshold:
        return "Not detected"
    num_detected = sum([e > low_expression_threshold for e in expr])
    if num_detected == 1:
        return "Detected in single"
    if num_detected < len(row) / 3:
        return "Detected in some"
    if num_detected < len(row):
        return "Detected in many"
    return "Detected in all"


# Adatiss specificity scores.


def adatiss(df, name_to_uberon_mapping):
    def _pack_adatiss_row(row, columns, name_to_uberon_mapping):
        """Given a row with Adatiss results for a given gene, pack them into the dictionary ready for output."""
        return [
            {
                "bodyPartLevel": "tissue",
                "bodyPartName": col,
                "bodyPartId": name_to_uberon_mapping[col],
                "adatissScore": round(val, 3),
            }
            for col, val in zip(columns, row.values.tolist())
        ]

    # Prepare Adatiss input file. Original column names cannot be ingested by Adatiss because they contain special
    # characters. Hence temporary identifiers are used for Adatiss.
    column_map = {original_name: f"tissue_{i}" for i, original_name in enumerate(df.columns)}
    reverse_column_map = {v: k for k, v in column_map.items()}
    df.rename(columns=column_map).to_csv("/tmp/adatiss_input.csv", index=True)

    # Adatiss also requires a sample to tissue name map file.
    with open("/tmp/adatiss_sample_phenotype.csv", "w") as outfile:
        outfile.write("sample_ID,tissue\n")
        for original_name, adatiss_name in column_map.items():
            outfile.write(f"{adatiss_name},{original_name}\n")

    # Calculate Adatiss specificity scores.
    subprocess.call(
        [
            "Rscript",
            os.path.join(os.path.dirname(__file__), "process_adatiss.R"),
            os.path.join(os.path.dirname(__file__), "AdaTiSS_fn.R"),
            "/tmp/adatiss_input.csv",
            "/tmp/adatiss_sample_phenotype.csv",
            "/tmp/adatiss_output.tsv",
        ]
    )

    # Read the results.
    adatiss_output = pd.read_table("/tmp/adatiss_output.tsv")
    adatiss_output.rename(columns=reverse_column_map, inplace=True)

    return adatiss_output.apply(
        functools.partial(
            _pack_adatiss_row, columns=adatiss_output.columns, name_to_uberon_mapping=name_to_uberon_mapping
        ),
        axis=1,
    )
