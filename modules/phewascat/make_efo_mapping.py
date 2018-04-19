'''
recreates mapping table

use --overwrite to force rewriting existing mappings
'''

import csv
import os
import logging
import argparse
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper

import click
import requests
from tqdm import tqdm
from ontoma import OnToma

from constants import *

# #ontoma's logger is useful to find out mapping issues
# from ontoma import logger as ontomalogger
# ontomalogger.setLevel(logging.INFO)

logger = logging.getLogger(__name__)

PHEWAS_CATALOG_URL = 'https://storage.googleapis.com/otar000-evidence_input/PheWAScatalog/phewas-catalog.csv'
PHEWAS_PHECODE_MAP_URL = 'https://phewascatalog.org/files/phecode_icd9_map_unrolled.csv.zip'

def download_ic9_phecode_map(url=PHEWAS_PHECODE_MAP_URL):
    '''phewascatalog maps ICD9 to a strange PheCode, but they provide mappings
    '''
    with requests.get(url) as phecode_res:
        # let us state clearly that I hate zip files. use gzip people!
        phecode_zip = ZipFile(BytesIO(phecode_res.content))
        # you can't just pass a ZipFile object as you would do a gZipfile
        # object. You actually have to extract it and to do that you need to
        # pass the name of the file it contains.
        phecode_file = phecode_zip.open(phecode_zip.namelist()[0])

    with TextIOWrapper(phecode_file) as phecode_map:
        reader = csv.DictReader(phecode_map)
        return {row['phecode']:row['icd9'] for row in reader}



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true',
                        help='force rewrite of existing mapping file')
    args = parser.parse_args()


    phecode_to_ic9 = download_ic9_phecode_map()
    otmap = OnToma()
    moduledir = os.path.dirname(os.path.realpath(__file__))
    efo_fn = os.path.join(moduledir, 'phewascat-diseaseterms.tsv')
    fieldnames = ['query', 'term', 'label', 'source', 'quality', 'action']

    if (not args.overwrite) and os.path.exists(efo_fn):
        logger.info('Reading existing mappings...')
        with open(efo_fn, 'r', newline='') as efo_f:
            eforeader = csv.DictReader(efo_f, fieldnames, delimiter='\t')
            #create a set of those terms that have been mapped succesfully before
            done = {row['query'] for row in eforeader if row['term']}

    efo_file = open(efo_fn, 'w', newline='') if args.overwrite else open(efo_fn, 'a', newline='')

    efowriter = csv.DictWriter(efo_file, fieldnames, delimiter='\t')
    if args.overwrite:
        logger.warning('Writing over existing files if any')
        efowriter.writeheader()

    #find EFO term
    mapped = 0
    processed = 0
    with requests.get(PHEWAS_CATALOG_URL, stream=True) as req:
        catalog = tqdm(csv.DictReader(req.iter_lines(decode_unicode=True)),
                        total=TOTAL_NUM_PHEWAS)
        for processed, row in enumerate(catalog):
            logger.debug(row)

            if done and (row['phewas_string'] in done):
                continue

            #first try using the ICD9 code that we are given
            try:
                efoid = otmap.find_term(phecode_to_ic9[row['phewas_code']],
                                        code="ICD9CM", verbose=True)
            except KeyError as e:
                logger.warning('Could not find %s in the phewascatalog phecode map', e)
                efoid = None

            if efoid:
                efoid['query'] = row['phewas_string']
                efowriter.writerow(efoid)
                mapped += 1
            else:
                # if we cannot find an ICD9 match, try to match for string using
                # OnToma
                efomatch = otmap.find_term(row['phewas_string'], verbose=True)
                if efomatch:
                    mapped += 1
                    efomatch['query'] = row['phewas_string']
                    efowriter.writerow(efomatch)
                else:
                    # no match gives you only the query and a blank line
                    efowriter.writerow({'query': row['phewas_string']})


    efo_file.close()
    click.echo("Completed. Parsed {} rows. "
               "Found {} EFOids. "
               "Skipped {} ".format(processed+1,
                                    mapped,
                                    processed+1-mapped))



if __name__ == '__main__':
    main()
