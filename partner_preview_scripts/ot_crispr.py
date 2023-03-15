#!/usr/bin/env python3
"""Parser to generate evidence for ot_crispr datasets.

This dataset consist of a series of OTAR projects studying various diseases using genome-wide crisp/cas9 knock-outs.
 - The results are expected to arrive in MAGeCK format.
 - The study level metadata is expected to come via filling out a Google spreadseet.
 - These spreadseet is then converted into a json and version-ed in the PPP-evidencie-configuration repository.
 - The generated evidence is exptected to be validated against the OT evidence schema.

"""
import argparse
from datetime import datetime
import gzip
import json
import logging
import sys

import pandas as pd

# The statisticalTestTail is inferred by the column name which is being filtered on:
FILTER_COLUMN_MAP = {
    "pos|fdr": "upper tail",
    "neg|fdr": "lower tail",
    "neg|p-value": "lower tail",
    "pos|p-value": "upper tail",
}


class OTAR_CRISPR_study_parser(object):
    def __init__(self, study_url: str) -> None:

        self.study_file = study_url

        # Store the study dataframe after dropping problematic studies:
        self.study_df = (
            pd.read_json(study_url)
            # drop rows with no study id or data file:
            .loc[lambda df: df.studyId.notna() & df.dataFiles.notna()]
        )

        # Test and warn if multiple studies have the same study id:
        duplicated_study_ids = self.study_df.loc[
            lambda df: df.studyId.duplicated()
        ].studyId.tolist()
        assert (
            len(duplicated_study_ids) == 0
        ), f'Multiple studies have the same study id: {", ".join(duplicated_study_ids)}'

        logging.info(
            f"Number of studies processed: {len(self.study_df.studyId.unique())}"
        )

        projects = self.study_df.projectId.unique()
        logging.info(f'Number of projects: {len(projects)} ({", ".join(projects)})')

    def generate_evidence(self, data_folder: str) -> None:

        # Looping through the studies and generating evidence:
        # Reading all data files and filter for significant hits:
        study_columns = [
            "releaseDate",
            "releaseVersion",
            "studyId",
            "dataFiles",
            "dataFileType",
            "filterColumn",
            "threshold",
            "projectId",
            "ControlDataset",
            "projectDescription",
        ]

        # hits is a pd.Series with pd.DataFrames as values.
        hits = (
            self.study_df[study_columns]
            .explode("dataFiles")
            .assign(
                dataFile=lambda df: df.apply(
                    lambda x: f'{data_folder}/{x["projectId"]}/{x["dataFiles"]}', axis=1
                )
            )
            .assign(
                ControlDataset=lambda df: df.apply(
                    lambda x: f'{data_folder}/{x["projectId"]}/{x["ControlDataset"]}'
                    if pd.notna(x["ControlDataset"])
                    else None,
                    axis=1,
                )
            )
            # TODO: parsing the data files should be file type dependent!
            # The following apply returns pd.DataFrames:
            .apply(self.parse_MAGeCK_file, axis=1)
        )

        # Concatenate all hits into one single dataframe:
        hits_df = (
            pd.concat(hits.to_list()).reset_index(drop=True)
            # Cleaning gene id column:
            .assign(
                targetFromSourceId=lambda df: df.targetFromSourceId.apply(
                    self.cleaning_gene_id
                )
            )
        )

        # Merging:
        evidence_fields = [
            "targetFromSourceId",
            "diseaseFromSourceMappedId",
            "projectDescription",
            "projectId",
            "studyId",
            "studyOverview",
            "contrast",
            "crisprScreenLibrary",
            "cellType",
            "cellLineBackground",
            "geneticBackground",
            "statisticalTestTail",
            "resourceScore",
            "log2FoldChangeValue",
            "releaseDate",
            "releaseVersion",
        ]
        self.merged_dataset = (
            self.study_df.assign(
                statisticalTestTail=lambda df: df.filterColumn.map(FILTER_COLUMN_MAP)
            )
            .merge(hits_df, on="studyId", how="inner")
            .explode("diseaseFromSourceMappedId")
            .filter(items=evidence_fields)
            .assign(datasourceId="ot_crispr", datatypeId="ot_partner")
        )

    @staticmethod
    def cleaning_gene_id(gene_id: str) -> str:
        """
        Expandable set of string processing steps to clean gene identifiers"""
        
        Examples:
            >>> cleaning_gene_id("ENSG00000187123_LYPD6")
            >>> "ENSG00000187123"
        """

        # ENSG00000187123_LYPD6 -> ENSG00000187123
        gene_id = gene_id.split("_")[0]

        return gene_id

    @staticmethod
    def parse_MAGeCK_file(row: pd.Series) -> pd.DataFrame:
        """This function returns a pandas dataframe with the datafile and with properly named columns"""

        datafile = row["dataFile"]
        filterColumn = row["filterColumn"]
        threshold = float(row["threshold"])
        studyId = row["studyId"]
        controlDataFile = row["ControlDataset"]

        # Read data, filter and rename columns:
        mageck_df = (
            pd.read_csv(datafile, sep="\t")
            .rename(
                columns={
                    filterColumn: "resourceScore",
                    "id": "targetFromSourceId",
                    "neg|lfc": "log2FoldChangeValue",
                }
            )
            .loc[lambda df: df.resourceScore <= threshold][
                ["targetFromSourceId", "resourceScore", "log2FoldChangeValue"]
            ]
            .assign(studyId=studyId)
        )

        # Applying control if present:
        if pd.isna(controlDataFile):
            logging.info(f"Number of genes reach threshold: {len(mageck_df)}")
            return mageck_df

        # Read control data, filter and rename columns:
        logging.info(f"Reading control data file: {controlDataFile}")
        controlHits = (
            pd.read_csv(controlDataFile, sep="\t")
            .rename(columns={filterColumn: "resourceScore", "id": "targetFromSourceId"})
            .loc[lambda df: df.resourceScore <= threshold]["targetFromSourceId"]
            .tolist()
        )

        # Excluding control genes:
        mageck_df = mageck_df.loc[lambda df: df.targetFromSourceId.isin(controlHits)]

        logging.info(f"Number of genes reach threshold: {len(mageck_df)}")
        return mageck_df

    def write_evidence(self, output_file: str) -> None:
        """Write the merged evidence to file"""

        json_list = [
            json.dumps(row.dropna().to_dict())
            for _, row in self.merged_dataset.iterrows()
        ]
        with gzip.open(output_file, "wt") as f:
            f.write("\n".join(json_list) + "\n")


def main(study_table, output_file, data_folder) -> None:
    # Parsing study table:
    parser = OTAR_CRISPR_study_parser(study_table)

    # Read data:
    parser.generate_evidence(data_folder)

    # Save data:
    parser.write_evidence(output_file)


if __name__ == "__main__":

    # Parsing arguments:
    parser = argparse.ArgumentParser(
        description="Script to parse internally generated datasets for OTAR projects."
    )

    parser.add_argument(
        "-s",
        "--study_table",
        type=str,
        required=True,
        help="A JSON file with study level metadata.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=str,
        help="Output json.gz file with the evidece.",
    )
    parser.add_argument(
        "-l",
        "--log_file",
        required=False,
        type=str,
        help="Logs are saved into this file.",
    )
    parser.add_argument(
        "-d",
        "--data_folder",
        required=True,
        type=str,
        help="Folder with the data files.",
    )
    args = parser.parse_args()

    study_table = args.study_table
    log_file = args.log_file
    data_folder = args.data_folder

    if not args.output:
        output_file = datetime.today().strftime("otar_crispr-%Y-%m-%d.json.gz")
    else:
        output_file = args.output

    # Configure logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Report input data:
    logging.info(f"Study information is read from: {study_table}")
    logging.info(f"Evidence saved to: {output_file}")
    logging.info(f"Data files read from: {data_folder}")

    # Process data:
    main(study_table, output_file, data_folder)
