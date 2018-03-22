import logging
import requests
import csv
import json
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper
import obonet
import python_jsonschema_objects as pjs 
from common.HGNCParser import GeneParser
from ontoma import OnToma
import click

#ontoma's logger is useful to find out mapping issues
from ontoma import logger as ontomalogger
ontomalogger.setLevel(logging.DEBUG)

#this module's logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


PHEWAS_CATALOG_URL = 'https://storage.googleapis.com/otar000-evidence_input/PheWAScatalog/phewas-catalog.csv'
PHEWAS_PHECODE_MAP_URL = 'https://phewascatalog.org/files/phecode_icd9_map_unrolled.csv.zip'




'''import the json schema and create python objects with a built-in validator
'''
click.echo('downloading schema')
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
click.echo('python objects created from schema')

'''define once where evidence is coming from
'''
prov = provenance(literature={"references":[
                    {"lit_id":"http://europepmc.org/articles/PMC3969265"}]}, 
                    database={"version":"2013-12-31T09:53:37+00:00", 
                    "id":"PHEWAS Catalog"}
                    )

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




'''gene symbol <=> ENSGID mappings'''

gene_parser = GeneParser()
logger.info('Parsing gene data from HGNC...')
gene_parser._get_hgnc_data_from_json()
ensgid = gene_parser.genes

'''prepare an object'''


PhewasEv = ot_objects['Genetics-basedEvidenceStrings']

click.echo('donwloading phecode mapping')
'''phewascatalog's Phecode <=> ICD9 mappings'''

phecode_to_ic9 = download_ic9_phecode_map()

click.echo('initializing ontoma')
otmap = OnToma()

skipped = []
built = 0
click.echo('begin processing phewascatalog evidences')
with open('output/phewas_test.json', 'w') as outfile:

    with requests.get(PHEWAS_CATALOG_URL, stream=True) as r:
        for i, row in enumerate(csv.DictReader(r.iter_lines(decode_unicode=True))):
            if i == 2:
                break
            logger.debug(row)
            pev = PhewasEv(type = 'genetic_association',
                        access_level = "public", 
                        sourceID = "phewas_catalog",
                        validated_against_schema_version = '1.2.8'
                        )
        
            '''find EFO term'''

            try:
                #first try using the ICD9 code that we are given
                efoid = otmap.find_efo(phecode_to_ic9[row['phewas_code']],
                                        code="ICD9CM")
            except KeyError as e:
                logger.error('No phecode <=> ICD9CM map for {}'.format(e))
                efoid = None

            if not efoid:
                logger.warning('Could not find a match '
                             'for {} in {} mappings. '.format(phecode_to_ic9[row['phewas_code']],"ICD9CM"))
                efoid = otmap.find_efo(row['phewas_string'])
            

            if efoid:
                logger.warning('Found {} for {}'.format(efoid,row['phewas_string']))
                pev['disease'] = make_disease(efoid)
            else: 
                logger.warning("Could not find disease: {}".format(row['phewas_string']))
                skipped.append((row['phewas_string'],phecode_to_ic9[row['phewas_code']]))
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

        skipped_uniq = set(skipped)
        logger.warning("Completed. Parsed {} rows. "
                       "Built {} evidences. " 
                       "Skipped {} ".format(i,built,len(skipped_uniq))
                       )

        with open('output/phewas_no_efo_codes.txt','wt',encoding='utf-8') as skipfile:
            for item in skipped_uniq:
                skipfile.write("\n".join(lines))

