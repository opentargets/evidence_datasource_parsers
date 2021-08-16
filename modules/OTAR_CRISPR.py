#!/usr/bin/env python3

import argparse
from datetime import datetime
import gzip
import json
import logging
import sys

import pandas as pd

FILTER_COLUMN_MAP = {
    'pos|fdr': 'upper tail',
    'neg|fdr': 'lower tail',
    'neg|p-value': 'lower tail',
    'pos|p-value': 'upper tail'
}

class OTAR_CRISPR_study_parser(object):
    def __init__(self, study_file: str) -> None:

        self.study_file = study_file

        # Read and process study table:
        study_df = pd.read_csv(study_file, sep='\t')

        # If present, the column description row is dropped:
        if 'Studied cell type' in study_df.loc[0].tolist():
            logging.info('Dropping column descriptions')
            study_df.drop(0, axis=0, inplace=True)

        logging.info(f'Number of studies processed: {len(study_df.studyId.unique())}')
        self.study_df = (
            study_df

            # drop rows with no study id or data file:
            .loc[study_df.studyId.notna() & study_df.dataFile.notna()]
            .assign(
                # Split mapped diseases to list:
                diseaseFromSourceMappedId=study_df.diseases.str.replace(' ', '').str.split('|'),

                # Split data files to list:
                dataFiles=study_df.dataFile.str.replace(' ', '').str.split('|')
            )

            # rename columns:
            .rename(
                columns={
                    'isDerived': 'isCellTypeDerived',
                    'library': 'crisprScreenLibrary',
                    'mode': 'crisprStudyMode',
                    'populations': 'contrast',
                    'studyDescription': 'studyOverview'
                }
            )
        )

    def generate_evidence(self, data_folder) -> None:
        """Processing the study dataframe"""

        # Looping through the studies and generating evidence:
        # Reading all data files and filter for significant hits:
        hits = (
            self.study_df[['studyId', 'dataFiles', 'dataFileType', 'filterColumn', 'threshold', 'projectId']]
            .explode('dataFiles')
            .assign(
                dataFile=lambda df: df.apply(lambda x: f'{data_folder}/{x["projectId"]}/{x["dataFiles"]}', axis=1)
            )
            # TODO: parsing the data files should be file type dependent:
            .apply(self.parse_MAGeCK_file, axis=1)
        )

        # Concatenate all hits into one single dataframe:
        hits_df = (
            pd.concat(hits.to_list())
            .reset_index(drop=True)
        )

        # Merging:
        self.merged_dataset = (
            self.study_df
            .assign(direction=lambda df: df.filterColumn.map(FILTER_COLUMN_MAP))
            .merge(hits_df, on='studyId', how='inner')
            .explode('diseaseFromSourceMappedId')
            .assign(
                datasourceId='ot_crispr',
                datatypeId='ot_partner'
            )
            [[
                'targetFromSourceId', 'diseaseFromSourceMappedId',
                'projectId', 'studyId', 'studyOverview', 'contrast', 'crisprScreenLibrary',
                'cellType', 'cellLineBackground', 'geneticBackground',
                'direction', 'resourceScore',
                'datasourceId', 'datatypeId'
            ]]
        )

    @staticmethod
    def parse_MAGeCK_file(row) -> pd.DataFrame:
        """This function returns with a pandas dataframe with the datafile and with properly named columns"""

        datafile = row['dataFile']
        filterColumn = row['filterColumn']
        threshold = float(row['threshold'])
        studyId = row['studyId']

        # Read data, filter and rename columns:
        try:
            mageck_df = (
                pd.read_csv(datafile, sep='\t')
                .rename(columns={filterColumn: 'resourceScore', 'id': 'targetFromSourceId'})
                .loc[lambda df: df.resourceScore <= threshold]
                [['targetFromSourceId', 'resourceScore']]
                .assign(studyId=studyId)
            )
            logging.info(f'Number of genes reach threshold: {len(mageck_df)}')
            return mageck_df
        except FileNotFoundError:
            logging.info(f'Study skipped as file was not found: {datafile}')

    def get_evidence_data(self) -> pd.DataFrame:
        """Returning the merged dataset"""
        return self.merged_dataset

    def write_evidence(self, output_file: str) -> None:
        """Write the merged evidence to file"""

        json_list = [json.dumps(row.dropna().to_dict()) for _, row in self.merged_dataset.iterrows()]
        with gzip.open(output_file, 'wt') as f:
            f.write('\n'.join(json_list) + '\n')


def main(study_table, output_file, data_folder) -> None:
    # Parsing study table:
    parser = OTAR_CRISPR_study_parser(study_table)

    # Read data:
    parser.generate_evidence(data_folder)

    # Save data:
    parser.write_evidence(output_file)


if __name__ == '__main__':

    # Parsing arguments:
    parser = argparse.ArgumentParser(description='Script to parse internally generated datasets for OTAR projects.')
    parser.add_argument('-s', '--study_table', help='A tsv file with the study description of the available projects.',
                        type=str, required=True)
    parser.add_argument('-o', '--output', help='Output json.gz file with the evidece.', required=False, type=str)
    parser.add_argument('-l', '--log_file', help='Generated logs can be saved into this file.',
                        required=False, type=str)
    parser.add_argument('-d', '--data_folder', help='Folder with the data files.', required=True, type=str)
    args = parser.parse_args()

    study_table = args.study_table
    log_file = args.log_file
    data_folder = args.data_folder

    if not args.output:
        output_file = datetime.today().strftime('otar_crispr-%Y-%m-%d.json.gz')
    else:
        output_file = args.output

    # Configure logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Report input data:
    logging.info(f'Study information is read from: {study_table}')
    logging.info(f'Evidence saved to: {output_file}')
    logging.info(f'Data files read from: {data_folder}')

    # Process data:
    main(study_table, output_file, data_folder)
