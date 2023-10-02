"""Ingests GTEx V8 data and computes the baseline expression data."""

import argparse
import functools
import gzip
import json

import pandas as pd

from . import metrics


class BaselineExpression:
    def read_gtex_data(self):
        """Read and preprocess GTEx source data."""
        # Read the source data. It's in the GCT format, which is a regular TSV, but the first two lines are metadata which we don't need.
        self.df = pd.read_csv(self.gtex_source_data, sep="\t", skiprows=2)
        # Drop column with gene names, which we don't need at all.
        self.df.drop(columns=["Description"], inplace=True)
        # Drop Ensembl gene version: "ENSG00000223972.2" → "ENSG00000223972".
        self.df["Name"] = df["Name"].str.split(".").str[0]
        # Ensembl gene IDs are used as index throughout the entire processing code.
        self.df.set_index(keys=["Name"], inplace=True)

    def get_name_to_uberon_mapping(self):
        """Fetch and process tissue name to UBERON ID lookup table."""
        # Load the tissue name to UBERON conversion file.
        df_name_mapping = pd.read_table(self.tissue_name_to_uberon_mapping_path)
        self.name_to_uberon_lut = {
            gtex_tissue: ontology_code
            for gtex_tissue, ontology_code in df_name_mapping[["gtex_tissue_name", "ontology_code"]].values.tolist()
        }
        # Verify that all columns are present in the mapping file.
        for tissue in self.df.columns:
            if tissue not in self.name_to_uberon_lut:
                raise AssertionError(f"GTEx tissue {tissue} is not present in the provided mappings file!")

    def filter_out_low_expression_genes(self):
        """Genes with very low expression tend to generate noisy data. One particular case is an expression vector where one tissue has expression of 0.1 TPM and all the other ones are 0.0 exactly. This causes the gene to be called as extremely specifically expressed, even though most of these issues are just false discoveries.
        A common threshold used for filtering out genes with very low expression is to filter out those where expression in all tissues is less than 1.0 TPM. This is also used by HPA, and supported by internal investigations into quality control.
        """
        max_expr_across_tissues = self.df.max(axis=1)
        expr_level_mask = max_expr_across_tissues >= self.low_expression_threshold
        self.df = self.df[expr_level_mask]

    def calculate_specificity_metrics(self):
        """Calculate and store several expression specificity metrics."""
        self.specificity = self.df.iloc[:, :0].copy()
        self.specificity["gini"] = self.df.apply(metrics.gini_coefficient, axis=1).round(3)
        self.specificity["hpaSpecificity"] = self.df.apply(
            functools.partial(metrics.hpa_specificity, low_expression_threshold=self.low_expression_threshold), axis=1
        )
        self.specificity["hpaDistribution"] = self.df.apply(
            functools.partial(metrics.hpa_distribution, low_expression_threshold=self.low_expression_threshold), axis=1
        )
        self.specificity["adatissScores"] = self.df.apply(
            functools.partial(metrics.adatiss, name_to_uberon_mapping=self.name_to_uberon_lut), axis=1
        )

    def pack_data_for_output(self):

        def _remove_adatiss_if_none(d):
            if d["adatissScores"] != d["adatissScores"]:  # using the fact that NaN != NaN
                return {k: v for k, v in d.items() if k != "adatissScores"}
            else:
                return d

        # Create a list to store JSON objects.
        json_list = []

        # Iterate over each row index.
        for idx in self.df.index:
            # Create a dictionary for the current row.
            row_dict = {
                "ensemblGeneId": idx,
                "expression": [
                    {
                        "bodyPartLevel": "tissue",
                        "bodyPartId": key,
                        "bodyPartName": self.name_to_uberon_lut[key],
                        "tpm": value,
                    }
                    for key, value in self.df.loc[idx].to_dict().items()
                ],
                "expressionSpecificity": _remove_adatiss_if_none(self.specificity.loc[idx].to_dict()),
            }
            # Append the dictionary to the list
            json_list.append(row_dict)

        # Output the data.
        with gzip.open(self.output_file_path, "wt", compresslevel=9) as f:
            for obj in json_list:
                f.write(json.dumps(obj))
                f.write("\n")

    def main(self):
        self.read_gtex_data()
        self.get_name_to_uberon_mapping()
        self.filter_out_low_expression_genes()
        self.calculate_specificity_metrics()
        self.pack_data_for_output()

    def __init__(self, gtex_source_data, tissue_name_to_uberon_mapping, output_file_path, low_expression_threshold):
        self.gtex_source_data = gtex_source_data
        self.tissue_name_to_uberon_mapping = tissue_name_to_uberon_mapping
        self.output_file_path = output_file_path
        self.low_expression_threshold = low_expression_threshold


parser = argparse.ArgumentParser()
parser.add_argument("--gtex-source-data", required=True, type=str, help="GTEx by-gene median TPM counts, in GCT format")
parser.add_argument(
    "--low-expression-threshold",
    required=True,
    type=float,
    help="Genes which are not expressed at least this much in at least one tissue will be filtered out",
)
parser.add_argument(
    "--tissue-name-to-uberon-mapping",
    required=True,
    type=str,
    help="A two column TSV file which maps GTEx tissue names to UBERON IDs.",
)
parser.add_argument(
    "--output-file-path",
    required=True,
    type=str,
    help="A path to output the output file. GZIP-compressed JSON, one object per line.",
)


if __name__ == "__main__":
    args = parser.parse_args()
    be = BaselineExpression(
        args.gtex_source_data, args.low_expression_threshold, args.tissue_name_to_uberon_mapping, args.output_file_path
    )
    be.main()
