'''tranform the phewascatalog.org CSV file in a set of JSON evidence objects
'''

import logging
import csv
import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np

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

def make_disease(disease_info, phecode):
    disease = {}
    disease["id"] = disease_info['efo_uri']
    disease['name'] = disease_info['efo_label']
    disease['source_name'] = disease_info['phewas_string'] + " [" + phecode + "]"
    return disease

def make_target(ensgid):
    target = {}
    target['target_type'] = "http://identifiers.org/cttv.target/gene_evidence"
    target['id'] = "http://identifiers.org/ensembl/{}".format(ensgid)
    target['activity'] = "http://identifiers.org/cttv.activity/unknown"
    return target


def make_variant(rsid, variant_id):
    variant = {}
    variant['type'] = "snp single"
    variant['rs_id'] = rsid
    if pd.notna(variant_id):
        variant["id"] = variant_id 
    return variant


def make_gene2variant(consequence_link):
    return {'provenance_type':PROVENANCE,
            'is_associated':True,
            'date_asserted':"2017-12-31T09:53:37+00:00",
            'evidence_codes':["http://purl.obolibrary.org/obo/ECO_0000205"],
            'functional_consequence': consequence_link
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

def ensgid_from_gene(gene, ensgid):
    ## find ENSGID ##
    gene = gene.strip("*")
    if gene not in ensgid:
        #__log__.debug("Skipping unmapped target %s",gene)
        #skipped_unmapped_target_cnt += 1
        return np.nan
    return ensgid[gene]


def write_variant_id(row, one2many_variants):
    if row["snp"] not in one2many_variants:
        variant_id = "{}_{}_{}_{}".format(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        return variant_id
    else:
        return np.nan

def main():
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

    ## load pheWAS data enriched from Genetics Portal

    phewas_w_consequences = pd.read_csv(Config.PHEWAS_CATALOG_W_CONSEQUENCES)
    phewas_w_consequences.rename(columns = {'rsid':'snp', 'gene_id': 'gene'}, inplace = True)

    ## gene symbol <=> ENSGID mappings ##
    gene_parser = GeneParser()
    __log__.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes


    built = 0
    skipped_phewas_string_cnt = 0
    skipped_unmapped_target_cnt = 0
    skipped_zero_length_cnt = 0
    __log__.info('Begin processing phewascatalog evidences..')
    with open(Config.PHEWAS_CATALOG_EVIDENCE_FILENAME, 'w') as outfile:

        with open(Config.PHEWAS_CATALOG_FILENAME) as r:
            #catalog = tqdm(csv.DictReader(r), total=TOTAL_NUM_PHEWAS)
            catalog = pd.read_csv(r, dtype={"chromosome":"str", "basepair":"str", "phewas_code": "str", "cases":"Int64", "odds_ratio":"float64", "p":"float64"})
            
            # Parsing genes
            catalog["gene"] = catalog["gene"].dropna().apply(lambda X: ensgid_from_gene(X, ensgid))
            catalog.dropna(subset=["gene"], inplace=True)
            # Merging dataframes: more records due to 1:many associations
            enriched_catalog = pd.merge(catalog, phewas_w_consequences, on=["gene", "snp"])

            one2many_variants = {i for i in phewas_w_consequences["snp"] if phewas_w_consequences["snp"].tolist().count(i) > 1}
            enriched_catalog["variant_id"] = enriched_catalog.dropna(how="any", subset=["chrom", "pos", "ref", "alt"]).apply(lambda X: write_variant_id(X, one2many_variants), axis=1)
            enriched_catalog.drop(["csq_arr", "most_severe_gene_csq"], axis=1, inplace=True)

            # Dropping exploded duplicates
            columns = ['chromosome', 'basepair', 'gene','snp', 'phewas_code', 'phewas_string', 'cases', 'odds_ratio','p','chrom','pos','most_severe_csq','consequence_link']
            enriched_catalog.drop_duplicates(subset=columns, inplace=True)

            # Only use data with p-value<0.
            enriched_catalog = enriched_catalog[enriched_catalog['p'] < 0.05]

            for i, row in enriched_catalog.iterrows():
                __log__.debug(row)

                gene = row["gene"]
                phewas_string = row['phewas_string'].strip()
                phewas_code = row['phewas_code'].strip()
                snp = row['snp']
                row_p = row['p']
                odds_ratio = row['odds_ratio']
                cases = row['cases']
                variant_id = row['variant_id']
                if pd.isna(row["consequence_link"]):
                    consequence_link = "http://purl.obolibrary.org/obo/SO_0001060"
                    functional_csq = "sequence_variant"
                else:
                    consequence_link = row["consequence_link"]
                    functional_csq = row["most_severe_csq"]
            

                pev = {}
                pev['type'] = 'genetic_association'
                pev['access_level'] = 'public'
                pev['sourceID'] = 'phewas_catalog'
                pev['validated_against_schema_version'] = '1.7.5'


                pev['target'] = make_target(gene)

                pev["variant"]= make_variant(snp, variant_id)
                pev["evidence"] = { "variant2disease": make_variant2disease(row_p, odds_ratio, cases),
                                    "gene2variant": make_gene2variant(consequence_link)}

                ## find EFO term ##
                if phewas_string not in mappings:
                    __log__.debug("Skipping unmapped disease %s", phewas_string)
                    skipped_phewas_string_cnt += 1
                    continue
                
                for disease in mappings[phewas_string]:
                    if disease['efo_uri'] != "NEW TERM REQUEST":
                        pev['disease'] = make_disease(disease, phewas_code)

                        # Evidence strings are unique based on the target, disease EFO term and Phewas string
                        pev['unique_association_fields'] = {'target_id': gene,
                                                            'disease_id' : pev['disease']['id'],
                                                            'phewas_string_and_code' : phewas_string + " [" + phewas_code + "]",
                                                            'variant_id': variant_id if pd.notna(variant_id) else snp}

                    
                        outfile.write("%s\n" % json.dumps(pev,
                            sort_keys=True, separators = (',', ':')))

                built +=1


        __log__.info("Completed. Parsed %s rows. "
            "Built %s evidences. "
            "Skipped %s ",i,built,i-built
            )
        __log__.debug("Skipped unmapped PheWAS string: {} \n Skipped unmapped targets: {} \n Skipped zero length gene name: {}".format(skipped_phewas_string_cnt, skipped_unmapped_target_cnt, skipped_zero_length_cnt))


if __name__ == '__main__':
    main()