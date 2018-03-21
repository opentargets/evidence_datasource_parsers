import requests
import csv
import json
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper
import obonet
from settings import Config
import python_jsonschema_objects as pjs 
from common.HGNCParser import GeneParser
from ontoma import OnToma


import logging
logger = logging.getLogger(__name__)

'''import the json schema and create python objects with a built-in validator
'''

schema = requests.get('https://raw.githubusercontent.com/opentargets/json_schema/ep-fixrelative/src/genetics.json')
logger.debug("Downloaded the following schema: {}".format(schema))

builder = pjs.ObjectBuilder(schema.json())
ot_objects = builder.build_classes()
logger.debug('The schema contains {}'.format([str(x) for x in dir(ot_objects)]))

target = ot_objects['Target']
disease = ot_objects['Disease']
variant = ot_objects['Variant<anonymous>']
gene2variant = ot_objects['Gene2variant<anonymous>']
variant2disease = ot_objects['Variant2disease<anonymous>']
provenance = ot_objects['ProvenanceType']

'''define once where evidence is coming from
'''
prov = provenance(literature={"references":[{"lit_id":"http://europepmc.org/articles/PMC3969265"}]}, 
                    database={"version":"2013-12-31T09:53:37+00:00", "id":"PHEWAS Catalog"})


def download_ic9_phecode_map(url=Config.PHEWAS_PHECODE_MAP_URL):
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
    return gene2variant(provenance_type=prov,
            is_associated=True, 
            date_asserted="2017-12-31T09:53:37+00:00",
            evidence_codes=["http://purl.obolibrary.org/obo/ECO_0000205"],
            functional_consequence='http://purl.obolibrary.org/obo/SO_0001060'
    )

def make_variant2disease(pval):
    return variant2disease(unique_experiment_reference='http://europepmc.org/articles/PMC3969265',
                           provenance_type=prov,
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



def main():

    '''gene symbol <=> ENSGID mappings'''

    gene_parser = GeneParser()
    logger.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes

    '''prepare an object'''

    
    PhewasEv = ot_objects['Genetics-basedEvidenceStrings']

    '''phewascatalog's Phecode <=> ICD9 mappings'''

    phecode_to_ic9 = download_ic9_phecode_map()

    otmap = OnToma()

    skipped = []
    built = 0
    with open('output/phewas_test.json', 'w') as outfile:

        with requests.get(Config.PHEWAS_CATALOG_URL, stream=True) as r:
            for i, row in enumerate(csv.DictReader(r.iter_lines(decode_unicode=True))):
                if i == 5000:
                    break
                logger.debug(row)
                pev = PhewasEv(type = 'genetic_association',
                            access_level = "public", 
                            sourceID = "phewas_catalog",
                            validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                            )
                try:
                    efoid = otmap.find_efo(phecode_to_ic9[row['phewas_code']],
                                            code="ICD9CM")
                except KeyError as e:
                    logger.error('No phecode map for {}'.format(e))
                    efoid = None
                    pass

                if not efoid:
                    try:
                        efoid = otmap.find_efo(row['phewas_string'])
                        pev['disease'] = make_disease(efoid)
                    except KeyError as e:
                        logger.warning('Could not find EFO ID for {}'.format(e))
                        skipped.append(e)
                        continue

                try:
                    pev['target'] = make_target(ensgid[row['gene'].strip('*')])
                except KeyError as e:
                    logger.warning("Could not find gene: {}".format(row['gene']))
                    skipped +=1
                    continue

                pev["variant"]= make_variant(row['snp'])
                pev["evidence"] = { "variant2disease": make_variant2disease(row['p']),
                                    "gene2variant": make_gene2variant() }
                
                pev['unique_association_fields'] = {'odds_ratio':row['odds_ratio'], 
                                                    'cases' : row['cases'], 
                                                    'phenotype' : row['phewas_string']}

                built +=1
                outfile.write("%s\n" % pev.serialize())

            logger.info("Completed. Parsed {} rows. Built {} evidences. Skipped {}".format(i,built,skipped))
            with open('output/phewas_no_efo_codes.txt') as skipfile:
                json.dump(skipped,skipfile)

    return

if __name__ == '__main__':
    main()
