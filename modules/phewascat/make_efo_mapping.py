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
from pathlib import Path

import click
import requests
from tqdm import tqdm
from constants import *

from ontoma import OnToma
#ontoma's logger is useful to find out mapping issues, but it's a bit loud
from ontoma import logger as ontomalogger
ontomalogger.setLevel(logging.WARNING)

__log__ = logging.getLogger(__name__)
__moduledir__ = Path(__file__).resolve().parent
__modulename__ = Path(__file__).parent.name


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
    otmap = OnToma(exclude='zooma')
    mappings = {}

    efo_fn = str(__moduledir__ / 'phewascat-diseaseterms.tsv')
    ## TODO this should be a constant I import from ontoma
    fieldnames = ['query', 'term', 'label', 'source', 'quality', 'action']

    #TODO: check for mappings on github, rather then here.

    ## check if we have mappings and if we should re-use them
    if (not args.overwrite) and os.path.exists(efo_fn):
        __log__.info('Caching existing mappings...')
        with open(efo_fn, 'r', newline='') as efo_f:
            eforeader = csv.DictReader(efo_f, fieldnames, delimiter='\t')
            #create a dict for those terms that have been already mapped
            # succesfully (or have been manually curated)
            mappings = {row['query']:row for row in eforeader if row['quality'] == 'match'}


    ## use local phewascatalog.org file if it exists where you run the script.
    if os.path.exists(os.path.join(str(__moduledir__), 'phewas-catalog.csv')):
        __log__.warning('Reading phewascatalog.org local file')
        catalog = tqdm(csv.DictReader(open(os.path.join(str(__moduledir__), 'phewas-catalog.csv'))),
                       total=TOTAL_NUM_PHEWAS)
    else:
        session = requests.Session()
        req = session.get(PHEWAS_CATALOG_URL, stream=True)
        catalog = tqdm(csv.DictReader(req.iter_lines(decode_unicode=True)),
                       total=TOTAL_NUM_PHEWAS)

    ## do some mapping now
    mappings = {}
    for processed, row in enumerate(catalog):
        __log__.debug(row)

        # skip duplicates disease. we need to find the mapping only once for all
        # evidences connected to that disease
        if row['phewas_string'] in mappings: continue

        #first try using the ICD9 code that we are given
        try:
            efoid = otmap.find_term(phecode_to_ic9[row['phewas_code']],
                                    code="ICD9CM", verbose=True)
        except KeyError as e:
            __log__.warning('Could not find %s in the phewascatalog phecode map', e)
            efoid = None

        if efoid:
            mappings[row['phewas_string']] = {**efoid, 'query': row['phewas_string']}
        else:
            # if we cannot find an ICD9 match, try to match for string using
            # OnToma
            efomatch = otmap.find_term(row['phewas_string'],
                                        verbose=True, suggest=True)
            if efomatch:
                mappings[row['phewas_string']] = {**efomatch, 'query': row['phewas_string']}
            else:
                # no match gives you only the query and a blank line
                mappings[row['phewas_string']] = {'query': row['phewas_string']}

    ## write or append mappings to a file
    efo_file = open(efo_fn, 'w', newline='') if args.overwrite else open(efo_fn, 'a', newline='')
    efowriter = csv.DictWriter(efo_file, fieldnames, delimiter='\t')
    if args.overwrite:
        __log__.warning('Writing over existing files if any')
        efowriter.writeheader()
    efowriter.writerows(mappings.values())
    efo_file.close()

    click.echo("Completed. Parsed {} rows. "
               "Found {} EFOids. "
               "Skipped {} ".format(processed+1,
                                    len(mappings),
                                    processed+1-len(mappings)))



if __name__ == '__main__':
    main()
