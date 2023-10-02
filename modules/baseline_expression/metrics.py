"""Collection of expression specificity calculation metrics."""

import subprocess

import numpy as np
import pandas as pd


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


def adatiss(df):
    # Prepare Adatiss input file. Original column names cannot be ingested by Adatiss because they contain special characters. Hence temporary identifiers are used for Adatiss.
    column_map = {original_name: f"tissue_{i}" for i, original_name in enumerate(df.columns)}
    reverse_column_map = {v: k for k, v in column_map.items()}
    df.rename(columns=column_map).to_csv("adatiss_input.csv", index=True)

    # Adatiss also requires a sample to tissue name map file.
    with open("adatiss_sample_phenotype.csv", "w") as outfile:
        outfile.write("sample_ID,tissue\n")
        for original_name, adatiss_name in column_map.items():
            outfile.write(f"{adatiss_name},{original_name}\n")

    # Calculate Adatiss specificity scores.
    subprocess.call(["Rscript", "./process.R"])

    # Read the results.
    adatiss = pd.read_table("adatiss_output.tsv")
    adatiss.rename(columns=reverse_column_map, inplace=True)
    return adatiss
