import requests
import csv
import json
from zipfile import ZipFile
from io import BytesIO, TextIOWrapper
import obonet
from settings import Config
import python_jsonschema_objects as pjs 
from common.HGNCParser import GeneParser


import logging
logger = logging.getLogger(__name__)


def download_ic9_phecode_map(url=Config.PHEWAS_PHECODE_MAP_URL):
    with requests.get(url) as phecode_res:
        # let us state clearly that I hate zip files. use gzip people!
        phecode_zip = ZipFile(BytesIO(phecode_res.content))    
        # you can't just pass a ZipFile object as you would do a gZipfile
        # object. You actually have to extract it and to do that you need to
        # pass the name of the file it contains..     
        phecode_file = phecode_zip.open(phecode_zip.namelist()[0])

    with TextIOWrapper(phecode_file) as phecode_map:
        reader = csv.DictReader(phecode_map)
        return {row['phecode']:row['icd9'] for row in reader}


def make_uri(ontology_id):
    ontology_code = ontology_id.replace(':',"_")
    if ontology_code.startswith('EFO'):
        return {'id': 'http://www.ebi.ac.uk/efo/'+ontology_code}
    elif ontology_code.startswith('HP') or ontology_code.startswith('MP') :
        return {'id': 'http://purl.obolibrary.org/obo/' + ontology_code}
    elif ontology_code.startswith('Orphanet') :
        return {'id': 'http://www.orpha.net/ORDO/' + ontology_code}
    else:
        logger.error(ontology_code)
        raise Exception


def make_target_entity(ensgid):
    return {"target_type": "http://identifiers.org/cttv.target/gene_evidence", 
            "id": "http://identifiers.org/ensembl/{}".format(ensgid), 
            "activity": "http://identifiers.org/cttv.activity/unknown"}


def make_variant_entity(rsid):
    return {"type": "snp single", 
            "id": "http://identifiers.org/dbsnp/{}".format(rsid)}


def gene2variant_entity():
    return {'provenance_type': {
                "database":{"version":"2017-12-31T09:53:37+00:00","id":"PHEWAS Catalog",}
                },
            'is_associated': True, 
            'date_asserted' : "2017-12-31T09:53:37+00:00",
            'evidence_codes':["http://purl.obolibrary.org/obo/ECO_0000205"],
            'functional_consequence':'http://purl.obolibrary.org/obo/SO_0001060'
            }

def make_variant2disease_entity(pval):
    return {'unique_experiment_reference':'http://europepmc.org/articles/PMC3969265',
            'provenance_type': {
                        "literature": {
                            "references":[{"lit_id":"http://europepmc.org/articles/PMC3969265"}]
                            },
                        "database": {
                            "version":"2013-12-31T09:53:37+00:00",
                            "id":"PHEWAS Catalog"
                            }
                        },
            'is_associated': True,
            'resource_score':{
                'type': 'pvalue', 
                'method': {
                    "description":"pvalue for the phenotype to snp association."
                    },
                "value":float(pval)
                },
            'date_asserted': "2013-12-31T09:53:37+00:00",
            'evidence_codes': ['http://identifiers.org/eco/PheWAS'],
            }


class OTOntoMapper(object):
    '''Open Targets ontology mapping cascade

    if you have an id, find a xref to EFO
    elif search exact match to name in EFO (or synonyms)
    elif search fuzzy match to name in EFO
    elif search in OLS
    elif search in Zooma High confidence set

    '''
    def __init__(self):

        '''Parse the ontology obo files for exact match lookup'''

        self.efo = obonet.read_obo(Config.EFO_URL)
        logger.info('EFO parsed. Size: {} nodes'.format(len(self.efo)))
        self.hp = obonet.read_obo(Config.HP_URL)
        logger.info('HP parsed. Size: {} nodes'.format(len(self.hp)))

        
        '''Create name mappings'''

        # id_to_name = {id_: data['name'] for id_, data in efo.nodes(data=True)}
        self.name_to_efo = {data['name']: id_ 
                            for id_, data in self.efo.nodes(data=True)}
        
        logger.debug('We can now lookup "asthma" and get: {}'.format(self.name_to_efo['asthma']))

        self.name_to_hp = {data['name']: id_ 
                           for id_, data in self.hp.nodes(data=True)}
        
        logger.debug('We can now lookup "Phenotypic abnormality" and get: {}'.format(self.name_to_hp['Phenotypic abnormality']))

        '''Download OXO xrefs'''
        payload={"ids":[],"inputSource":"ICD9CM","mappingTarget":["EFO"],"distance":"3"}
         # HTTP POST https://www.ebi.ac.uk/spot/oxo/api/search?size=1000

    def hp_lookup(self, name):
        return self.name_to_hp[name]
    
    def efo_lookup(self, name):
        return self.name_to_efo[name]

    def oxo_lookup(self, other_ontology_id):
        '''should return an EFO code for any given xref'''
        return None


def main():

    '''gene symbol <=> ENSGID mappings'''

    gene_parser = GeneParser()
    logger.info('Parsing gene data from HGNC...')
    gene_parser._get_hgnc_data_from_json()
    ensgid = gene_parser.genes
    logger.debug(ensgid['BRAF'])


    '''prepare an object'''

    schema = requests.get('https://raw.githubusercontent.com/opentargets/json_schema/ep-fixrelative/src/genetics.json')

    logger.debug("Downloaded the following schema:")
    logger.debug(schema)

    builder = pjs.ObjectBuilder(schema.json())
    ns = builder.build_classes()
    logger.debug([str(x) for x in dir(ns)])
    PhewasEv = ns['Genetics-basedEvidenceStrings']
    logger.debug(PhewasEv(type='genetic_association'))

    '''Phecode <=> ICD9 mappings'''

    phecode_to_ic9 = download_ic9_phecode_map()
    logger.debug('Checking that the ic9 code for the 588 phecode is: %s' % phecode_to_ic9['588'])

    otmap = OTOntoMapper()
    logger.debug(otmap.efo_lookup('asthma'))


    with requests.get(Config.PHEWAS_CATALOG_URL, stream=True) as r:
        for i, row in enumerate(csv.DictReader(r.iter_lines(decode_unicode=True))):
            if i == 10:
                break
            logger.debug(row)
            pev = PhewasEv(type = 'genetic_association',
                           access_level = "public", 
                           sourceID = "phewas_catalog",
                           validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                           )
            
            efoid = otmap.oxo_lookup(phecode_to_ic9[row['phewas_code']])

            if not efoid:
                try:
                    efoid = otmap.efo_lookup(row['phewas_string'])
                    pev['disease'] = make_uri(efoid)
                except KeyError as e:
                    logger.warning('Could not find EFO ID for {}'.format(e))
                    continue

            logger.debug(pev)

            try:
                pev['target'] = make_target_entity(ensgid[row['gene'].strip('*')])
            except KeyError as e:
                logger.error(row['gene'])
                continue

            pev["variant"]= make_variant_entity(row['snp'])
            pev["evidence"] = { "variant2disease": make_variant2disease_entity(row['p']),
                                "gene2variant": gene2variant_entity() }
            
            pev['unique_association_fields'] = {'odds_ratio':row['odds_ratio'], 
                                                'cases' : row['cases'], 
                                                'phenotype' : row['phewas_string']}

    
            logger.debug(pev)


    return

if __name__ == '__main__':
    main()
