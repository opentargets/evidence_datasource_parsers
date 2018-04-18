'''tranform the phewascatalog.org CSV file in a set of JSON evidence objects
'''

import logging
import csv
import os
import sys
import json
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper

import requests
import obonet
from tqdm import tqdm
import python_jsonschema_objects as pjs
from python_jsonschema_objects.validators import ValidationError

from common.HGNCParser import GeneParser
from constants import *

#this module's logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


'''import the json schema and create python objects with a built-in validator
'''
logger.info('downloading schema')
schema = requests.get(PHEWAS_SCHEMA)
logger.debug("Downloaded the following schema: {}".format(schema))

builder = pjs.ObjectBuilder(schema.json())
ot_objects = builder.build_classes()
logger.debug('The schema contains {}'.format([str(x) for x in dir(ot_objects)]))

target = ot_objects['Target']
disease = ot_objects['Disease']
variant = ot_objects['Variant<anonymous>']
gene2variant = ot_objects['Gene2variant<anonymous>']
variant2disease = ot_objects['Variant2disease<anonymous>']
logger.info('python objects created from schema')




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
    return target(target_type = "http://identifiers.org/cttv.target/gene_evidence",
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

    ## load prepared mappings
    mappings = {}
    try:
        moddir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(moddir, 'terms.tsv')) as mapf:
            reader = csv.DictReader(mapf,delimiter='\t')
            mappings = {row['query']:row['term'] for row in reader}
    except FileNotFoundError:
        sys.exit('ABORTING. No existing mappings found. Please run make_efo_mapping.py first.')

    ## gene symbol <=> ENSGID mappings ##
    gene_parser = GeneParser()
    logger.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes

    ## phewascatalog's Phecode <=> ICD9 mappings ##
    phecode_to_ic9 = download_ic9_phecode_map()


    ## prepare an object ##
    PhewasEv = ot_objects['Genetics-basedEvidenceStrings']

    built = 0
    logger.info('Begin processing phewascatalog evidences..')
    with open(os.path.join(outputdir, 'phewas_catalog.json'), 'w') as outfile:

        with requests.get(PHEWAS_CATALOG_URL, stream=True) as r:
            for i, row in enumerate(tqdm(csv.DictReader(r.iter_lines(decode_unicode=True)),total=215107)):
                logger.debug(row)
                pev = PhewasEv(type = 'genetic_association',
                               access_level = "public",
                               sourceID = "phewas_catalog",
                               validated_against_schema_version = '1.2.8'
                              )

                '''find EFO term'''

                try:
                    pev['disease'] = make_disease(mappings[row['phewas_string']])
                except KeyError as e:
                    logger.error('No mapping for {}. Skipping evidence'.format(e))
                    continue
                except ValidationError:
                    logger.error('Empty mapping for {}. '
                                 'Skipping evidence'.format(row['phewas_string']))
                    continue


                ''' find ENSGID '''

                try:
                    pev['target'] = make_target(ensgid[row['gene'].strip('*')])
                except KeyError as e:
                    logger.error("Could not find gene: {}".format(row['gene']))
                    continue

                pev["variant"]= make_variant(row['snp'])
                pev["evidence"] = { "variant2disease": make_variant2disease(row['p']),
                                    "gene2variant": make_gene2variant() }

                pev['unique_association_fields'] = {'odds_ratio':row['odds_ratio'],
                                                    'cases' : row['cases'],
                                                    'phenotype' : row['phewas_string']}

                built +=1
                outfile.write("%s\n" % pev.serialize())

        logger.info("Completed. Parsed {} rows. "
                    "Built {} evidences. "
                    "Skipped {} ".format(i,built,i-built)
                    )


if __name__ == '__main__':
    main(outputdir=os.path.dirname(os.path.realpath(__file__)))