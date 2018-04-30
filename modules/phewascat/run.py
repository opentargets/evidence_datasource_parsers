'''tranform the phewascatalog.org CSV file in a set of JSON evidence objects
'''

import logging
import csv
import os
import sys
import json
import argparse
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper
from pathlib import Path


import requests
import obonet
from tqdm import tqdm
import python_jsonschema_objects as pjs
from python_jsonschema_objects.validators import ValidationError

from common.HGNCParser import GeneParser
from common.Utils import mapping_on_github, ghmappings
from constants import *


__log__ = logging.getLogger(__name__)
__moduledir__ = Path(__file__).resolve().parent
__modulename__ = Path(__file__).parent.name

#configure this module's logger
__log__.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
__log__.addHandler(ch)



'''import the json schema and create python objects with a built-in validator
'''
__log__.info('downloading schema')
schema = requests.get(PHEWAS_SCHEMA)
__log__.debug("Downloaded the following schema: %s", schema)

builder = pjs.ObjectBuilder(schema.json())
ot_objects = builder.build_classes()
__log__.debug('The schema contains %s', [str(x) for x in dir(ot_objects)])

target = ot_objects['Target']
disease = ot_objects['Disease']
variant = ot_objects['Variant<anonymous>']
gene2variant = ot_objects['Gene2variant<anonymous>']
variant2disease = ot_objects['Variant2disease<anonymous>']
__log__.info('python objects created from schema')


def download_ic9_phecode_map(url=PHEWAS_PHECODE_MAP_URL):
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


def make_disease(ontology_uri):
    return disease(id=ontology_uri)

def make_target(ensgid):
    return target(target_type="http://identifiers.org/cttv.target/gene_evidence",
                  id="http://identifiers.org/ensembl/{}".format(ensgid),
                  activity="http://identifiers.org/cttv.activity/unknown")


def make_variant(rsid):
    return variant(type="snp single",
                   id="http://identifiers.org/dbsnp/{}".format(rsid))


def make_gene2variant():
    return {'provenance_type':PROVENANCE,
            'is_associated':True,
            'date_asserted':"2017-12-31T09:53:37+00:00",
            'evidence_codes':["http://purl.obolibrary.org/obo/ECO_0000205"],
            'functional_consequence':'http://purl.obolibrary.org/obo/SO_0001060'
            }

def make_variant2disease(pval):
    return variant2disease(unique_experiment_reference='http://europepmc.org/articles/PMC3969265',
                           provenance_type=PROVENANCE,
                           is_associated=True,
                           date_asserted="2013-12-31T09:53:37+00:00",
                           evidence_codes=['http://identifiers.org/eco/PheWAS'],
                           resource_score={'type': 'pvalue',
                                           'method': {
                                           "description":"pvalue for the phenotype to snp association."
                                          },
                            "value":float(pval)
                            }
                           )



def main(outputdir):
    parser = argparse.ArgumentParser()
    parser.add_argument('--local', action='store_true',
                        help='force using local mapping instead of the github version')
    args = parser.parse_args()

    ## load prepared mappings
    mappings = {}
    if (not args.local) and mapping_on_github(__modulename__):
        __log__.info('Downloading prepared mappings from github')
        __log__.info('Using only mappings marked as `match`')
        with requests.get(ghmappings(__modulename__), stream=True) as rmap:
            for row in csv.DictReader(rmap.iter_lines(decode_unicode=True),delimiter='\t'):
                if row['quality'] and row['quality'].strip() == 'match':
                    mappings[row['query'].strip()] = row['term'].strip()
        __log__.info('Parsed %s mappings', len(mappings))

    try:
        moddir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(moddir, 'phewascat-diseaseterms.tsv')) as mapf:
            reader = csv.DictReader(mapf,delimiter='\t')
            mappings = {row['query']:row['term'] for row in reader
                        if (row['quality'] and row['quality'] == 'match')}
    except FileNotFoundError:
        sys.exit('ABORTING. No existing mappings found. Please run make_efo_mapping.py first.')

    ## gene symbol <=> ENSGID mappings ##
    gene_parser = GeneParser()
    __log__.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes

    ## phewascatalog's Phecode <=> ICD9 mappings ##
    phecode_to_ic9 = download_ic9_phecode_map()


    ## prepare an object ##
    PhewasEv = ot_objects['Genetics-basedEvidenceStrings']

    built = 0
    __log__.info('Begin processing phewascatalog evidences..')
    with open(os.path.join(str(outputdir),'phewas-catalog.evidenceobjs.json'), 'w') as outfile:

        with requests.get(PHEWAS_CATALOG_URL, stream=True) as r:
            catalog = tqdm(csv.DictReader(r.iter_lines(decode_unicode=True)),
                           total=TOTAL_NUM_PHEWAS)
            for i, row in enumerate(catalog):
                __log__.debug(row)
                pev = PhewasEv(type = 'genetic_association',
                               access_level = "public",
                               sourceID = "phewas_catalog",
                               validated_against_schema_version = '1.2.8'
                              )

                ## find EFO term ##

                try:
                    pev['disease'] = make_disease(mappings[row['phewas_string']])
                except KeyError as e:
                    __log__.error('No mapping for %s. Skipping evidence',e)
                    continue
                except ValidationError:
                    __log__.error('Empty mapping for %s. '
                                 'Skipping evidence',row['phewas_string'])
                    continue


                ## find ENSGID ##

                try:
                    pev['target'] = make_target(ensgid[row['gene'].strip('*')])
                except ValidationError:
                    __log__.error("Invalid gene: "
                                  "http://identifiers.org/ensembl/%s",
                                  ensgid[row['gene'].strip('*')])
                    continue
                except KeyError as e:
                    __log__.error("Could not find gene: %s",row['gene'])
                    #TODO: deal with the `intergenic` case by looking up the
                    # postgap table or the snp2gene table in Open Targets
                    continue

                pev["variant"]= make_variant(row['snp'])
                pev["evidence"] = { "variant2disease": make_variant2disease(row['p']),
                                    "gene2variant": make_gene2variant() }

                pev['unique_association_fields'] = {'odds_ratio':row['odds_ratio'],
                                                    'cases' : row['cases'],
                                                    'phenotype' : row['phewas_string']}

                built +=1
                outfile.write("%s\n" % pev.serialize())

        __log__.info("Completed. Parsed %s rows. "
                    "Built %s evidences. "
                    "Skipped %s ",i,built,i-built
                    )


if __name__ == '__main__':
    p = Path(__file__).parents
    main(outputdir=p[2] / 'output')