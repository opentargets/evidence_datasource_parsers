import pandas as pd
import requests
from bs4 import BeautifulSoup, UnicodeDammit
import argparse
import logging
import logging.config
import sys


'''
This script retrieves TEP (target enabling package) from the structural genomics consoritum.
'''

def retrieve_tep_table() -> pd.DataFrame:

    def get_gene_name(row):
        return row.findAll('td')[0].text.split('/')

    def get_description(row):
        return row.findAll('td')[1].text.strip()

    def get_url(row):
        gene_cell = row.findAll('td')[0]
        tep_url = gene_cell.find('a').get('href')

        if not tep_url.startswith('http'):
            tep_url = 'https://www.thesgc.org' + tep_url

        return tep_url

    def get_therapeutic_area(row):
        therapeutic_area = row.findAll('td')[2].text
        return therapeutic_area

    url = 'https://www.thesgc.org/tep'

    response = requests.get(url)

    html = response.text
    uhtml = UnicodeDammit(html)

    soup = BeautifulSoup(uhtml.unicode_markup, features='html.parser')

    TEP_table = soup.findAll('table')[-1]

    TEP_raw_data = []

    for row in TEP_table.find('tbody').findAll('tr'):
        TEP_raw_data.append({
            'targetFromSourceId': get_gene_name(row),
            'url': get_url(row),
            'therapeuticArea': get_therapeutic_area(row),
            'description': get_description(row)
        })

    return pd.DataFrame(TEP_raw_data).explode('targetFromSourceId')


def main(outputFile: str) -> None:

    # The TEP list compiled as dataframe:
    tep_list = retrieve_tep_table()
    logging.info(f'Number of TEPs retrieved: {len(tep_list)}')
    logging.info(f'Number of unique gene symbols in the tep: {len(tep_list.targetFromSourceId.unique())}')

    if tep_list.targetFromSourceId.duplicated().any():
        duplicated_targets = tep_list.loc[tep_list.targetFromSourceId.duplicated()]
        logging.error(f'The following symbols were not unique: {duplicated_targets}')
        raise ValueError('The target symbol list is expected to be unique! Check logs for deatils.')

    # Saving data:
    logging.info(f'Saving data to {outputFile}.')
    tep_list.to_json(outputFile, compression='infer', orient='records', lines=True)


if __name__ == '__main__':

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(description='This script fetches TEP data from Structural Genomics Consortium.')
    parser.add_argument('--output_file', '-o', type=str, help='Output file. gzipped JSON', required=True)
    parser.add_argument('--log_file', type=str, help='File into which the logs are saved', required=False)
    args = parser.parse_args()

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    main(args.output_file)
