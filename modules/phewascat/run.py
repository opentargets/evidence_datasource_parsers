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
from tqdm import tqdm

from common.HGNCParser import GeneParser
from common.Utils import mapping_on_github, ghmappings, DuplicateFilter

PHEWAS_CATALOG_URL = 'https://storage.googleapis.com/otar000-evidence_input/PheWAScatalog/phewas-catalog.csv'

'''define once where evidence is coming from
'''
PROVENANCE = {'literature': {
            "references":[{"lit_id":"http://europepmc.org/articles/PMC3969265"}]
            },
             "database":{
                 "version":"2013-12-31T09:53:37+00:00",
                "id":"PHEWAS Catalog"
                }
}

TOTAL_NUM_PHEWAS = 4269549


__log__ = logging.getLogger(__name__)
__moduledir__ = Path(__file__).resolve().parent
__modulename__ = Path(__file__).parent.name

#configure this module's logger
__log__.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
__log__.addHandler(ch)
dup_filter = DuplicateFilter()
__log__.addFilter(dup_filter)


def make_disease(ontology_uri):
    disease = {}
    disease["id"] = ontology_uri
    return disease

def make_target(ensgid):
    target = {}
    target['target_type'] = "http://identifiers.org/cttv.target/gene_evidence"
    target['id'] = "http://identifiers.org/ensembl/{}".format(ensgid)
    target['activity'] = "http://identifiers.org/cttv.activity/unknown"
    return target


def make_variant(rsid):
    variant = {}
    variant['type'] = "snp single"
    variant['id'] = "http://identifiers.org/dbsnp/{}".format(rsid)
    return variant


def make_gene2variant():
    return {'provenance_type':PROVENANCE,
            'is_associated':True,
            'date_asserted':"2017-12-31T09:53:37+00:00",
            'evidence_codes':["http://purl.obolibrary.org/obo/ECO_0000205"],
            'functional_consequence':'http://purl.obolibrary.org/obo/SO_0001060'
            }

def make_variant2disease(pval):
    variant2disease = {}
    variant2disease['unique_experiment_reference']='http://europepmc.org/articles/PMC3969265'
    variant2disease['provenance_type']=PROVENANCE
    variant2disease['is_associated']=True
    variant2disease['date_asserted']="2013-12-31T09:53:37+00:00"
    variant2disease['evidence_codes']=['http://identifiers.org/eco/PheWAS']
    variant2disease['resource_score']={
        'type': 'pvalue',
        'method': {
            "description":"pvalue for the phenotype to snp association."
            },
        "value":float(pval)
        }

    return variant2disease



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
    else:
        moddir = os.path.dirname(os.path.realpath(__file__))
        inputfile = os.path.join(moddir, 'phewascat-diseaseterms.tsv')
        if os.path.exists(inputfile) and os.path.isfile(inputfile):
            with open(inputfile) as mapf:
                reader = csv.DictReader(mapf,delimiter='\t')
                mappings = {row['query']:row['term'] for row in reader if (row['quality'] and row['quality'] == 'match')}
        else:
            raise RuntimeError("Unable to read {}".format(inputfile))

    ## gene symbol <=> ENSGID mappings ##
    gene_parser = GeneParser()
    __log__.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes


    built = 0
    __log__.info('Begin processing phewascatalog evidences..')
    with open(os.path.join(str(outputdir),'phewas-catalog.evidenceobjs.json'), 'w') as outfile:

        with requests.get(PHEWAS_CATALOG_URL, stream=True) as r:
            catalog = tqdm(csv.DictReader(r.iter_lines(decode_unicode=True)),
                           total=TOTAL_NUM_PHEWAS)
            for i, row in enumerate(catalog):
                __log__.debug(row)

                phewas_string = row['phewas_string']
                gene = row['gene'].strip('*')
                snp = row['snp']
                row_p = row['p']
                odds_ratio = row['odds_ratio']
                cases = row['cases']

                pev = {}
                pev['type'] = 'genetic_association'
                pev['access_level'] = 'public'
                pev['sourceID'] = 'phewas_catalog'
                pev['validated_against_schema_version'] = '1.2.8'

                ## find EFO term ##
                if phewas_string not in mappings:
                    __log__.debug("Skipping unmapped disease %s",phewas_string)
                    continue
                pev['disease'] = make_disease(mappings[phewas_string])


                ## find ENSGID ##
                if gene not in ensgid:
                    __log__.debug("Skipping unmapped target %s",gene)
                    continue
                if not len(gene) or not len(ensgid[gene]):
                    __log__.debug("Skipping zero length gene name")
                    continue
                pev['target'] = make_target(ensgid[gene])

                pev["variant"]= make_variant(snp)
                pev["evidence"] = { "variant2disease": make_variant2disease(row_p),
                                    "gene2variant": make_gene2variant() }

                pev['unique_association_fields'] = {'odds_ratio':odds_ratio,
                                                    'cases' : cases,
                                                    'phenotype' : phewas_string}

                built +=1
                outfile.write("%s\n" % json.dumps(pev, 
                    sort_keys=True, separators = (',', ':')))

        __log__.info("Completed. Parsed %s rows. "
            "Built %s evidences. "
            "Skipped %s ",i,built,i-built
            )


if __name__ == '__main__':
    p = Path(__file__).parents
    main(outputdir=p[2] / 'output')