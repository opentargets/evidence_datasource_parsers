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

from settings import Config, file_or_resource

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

def make_disease(disease_info):
    disease = {}
    disease["id"] = disease_info['efo_uri']
    disease['name'] = disease_info['efo_label']
    disease['source_name'] = disease_info['phewas_string']
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

def make_variant2disease(pval, odds_ratio, cases):
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

    variant2disease['odds_ratio'] = odds_ratio
    variant2disease['cases'] = cases
    return variant2disease



def main(outputdir):

    ## load prepared mappings
    mappings = {}
    if mapping_on_github('phewascat'):
        __log__.info('Downloading prepared mappings from github')
        with requests.get(ghmappings('phewascat'), stream=True) as rmap:
            for row in csv.DictReader(rmap.iter_lines(decode_unicode=True),delimiter='\t'):
                phewas_str = row['Phewas_string'].strip()
                if phewas_str not in mappings:
                    mappings[phewas_str] = []
                mappings[phewas_str].append({'efo_uri':row['EFO_id'].strip(), 'efo_label':row['EFO_label'].strip(), 'phewas_string':phewas_str})
        __log__.info('Parsed %s mappings', len(mappings))
    else:
        __log__.error('Trait mapping file not found on GitHub: https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv')
        sys.exit()

    ## gene symbol <=> ENSGID mappings ##
    gene_parser = GeneParser()
    __log__.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes


    built = 0
    skipped_phewas_string_cnt = 0
    skipped_unmapped_target_cnt = 0
    skipped_zero_length_cnt = 0
    skipped_non_significant_cnt = 0
    __log__.info('Begin processing phewascatalog evidences..')
    with open(os.path.join(str(outputdir),'phewas-catalog.evidenceobjs.json'), 'w+') as outfile:

        with open(Config.PHEWAS_CATALOG_FILENAME) as r:
            #catalog = tqdm(csv.DictReader(r.iter_lines(decode_unicode=True)),
            #               total=TOTAL_NUM_PHEWAS)
            catalog = tqdm(csv.DictReader(r), total=TOTAL_NUM_PHEWAS)
            for i, row in enumerate(catalog):
                __log__.debug(row)

                # Only use data with p-value<0.05
                if float(row['p']) < 0.05:

                    phewas_string = row['phewas_string'].strip()
                    gene = row['gene'].strip('*')
                    snp = row['snp']
                    row_p = row['p']
                    odds_ratio = row['odds_ratio']
                    cases = row['cases']

                    pev = {}
                    pev['type'] = 'genetic_association'
                    pev['access_level'] = 'public'
                    pev['sourceID'] = 'phewas_catalog'
                    pev['validated_against_schema_version'] = '1.6.7'

                    ## find ENSGID ##
                    if gene not in ensgid:
                        __log__.debug("Skipping unmapped target %s",gene)
                        skipped_unmapped_target_cnt += 1
                        #__log__.info("Unmapped target: {}".format(row))
                        continue
                    if not len(gene) or not len(ensgid[gene]):
                        __log__.debug("Skipping zero length gene name")
                        skipped_zero_length_cnt += 1
                        continue
                    pev['target'] = make_target(ensgid[gene])

                    pev["variant"]= make_variant(snp)
                    pev["evidence"] = { "variant2disease": make_variant2disease(row_p, odds_ratio, cases),
                                        "gene2variant": make_gene2variant()}

                    ## find EFO term ##
                    if phewas_string not in mappings:
                        __log__.debug("Skipping unmapped disease %s", phewas_string)
                        __log__.info("Skipped disease: {}".format(phewas_string))
                        # __log__.info("Type: {}".format(type(phewas_string)))
                        # __log__.info("Dict: {}".format(mappings['Intestinal infection']))
                        skipped_phewas_string_cnt += 1
                        continue

                    for disease in mappings[phewas_string]:
                        if disease['efo_uri'] != "NEW TERM REQUEST":
                            pev['disease'] = make_disease(disease)

                            # Evidence strings are unique based on the target, disease EFO term and Phewas string
                            pev['unique_association_fields'] = {'target_id': ensgid[gene],
                                                                'disease_id' : pev['disease']['id'],
                                                                'phewas_string' : phewas_string}

                            #outfile.write("%s\n" % json.dumps(pev,
                            #    sort_keys=True, separators = (',', ':')))

                    built +=1

                else:
                    skipped_non_significant_cnt += 1

        __log__.info("Completed. Parsed %s rows. "
            "Built %s evidences. "
            "Skipped %s ",i,built,i-built
            )
        __log__.info("Skipped unmapped PheWAS string: {} \n Skipped unmapped targets: {} \n Skipped zero length gene name: {} \n Skipped non-significant association: {}".format(skipped_phewas_string_cnt, skipped_unmapped_target_cnt, skipped_zero_length_cnt, skipped_non_significant_cnt))


if __name__ == '__main__':
    p = Path(__file__).parents
    main(outputdir=p[1] / 'output')