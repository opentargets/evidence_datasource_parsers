'''general settings that all parsers can share
'''

import os
from pathlib import Path
import pkg_resources as res

# from envparse import env, ConfigurationError

# def read_option(option, cast=None,
#                 **kwargs):

#     try:
#         default_value = kwargs.pop('default')
#     except KeyError:
#         default_value = None

#     try:
#         # reading the environment variable with envparse
#         return env(option, cast=cast, **kwargs)
#     except ConfigurationError:
#        return default_value


def file_or_resource(fname=None):
    '''get filename and check if in getcwd then get from
    the package resources folder
    '''
    filename = os.path.expanduser(fname)

    resource_package = __name__
    resource_path = '/'.join(('resources', filename))

    if filename is not None:
        abs_filename = os.path.join(os.path.abspath(os.getcwd()), filename) \
                       if not os.path.isabs(filename) else filename

        return abs_filename if os.path.isfile(abs_filename) \
            else res.resource_filename(resource_package, resource_path)

class Config:
    '''shared settings
    '''
    HOME_DIR = str(Path.home())

    # schema version
    VALIDATED_AGAINST_SCHEMA_VERSION = '1.2.8'

    #Ontologies
    EFO_URL = 'https://github.com/EBISPOT/efo/raw/v2018-01-15/efo.obo'
    HP_URL = 'http://purl.obolibrary.org/obo/hp.obo'

    # HGNC
    GENES_HGNC = 'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'

    # PROGENY
    PROGENY_FILENAME = file_or_resource('progeny_normalVStumor_opentargets.txt')
    PROGENY_EVIDENCE_FILENAME = HOME_DIR + '/progeny-20-05-2018.json'

    # UKBIOBANK
    UKBIOBANK_FILENAME = file_or_resource('ukbiobank.txt')
    UKBIOBANK_EVIDENCE_FILENAME = HOME_DIR + '/ukbiobank-30-04-2018.json'

    # SLAPEnrich
    SLAPENRICH_FILENAME = file_or_resource('slapenrich_opentargets.tsv')
    SLAPENRICH_EVIDENCE_FILENAME = HOME_DIR + '/otar001_slapenrich-19-12-2017.json'

    # Gene2Phenotype
    G2P_FILENAME = 'DDG2P.csv.gz'
    G2P_EVIDENCE_FILENAME = 'gene2phenotype.json'

    # Genomics England
    GE_PANEL_MAPPING_FILENAME = file_or_resource('genomicsenglandpanelapp_panelmapping.csv')
    GE_EVIDENCE_FILENAME = 'genomics_england.json'
    GE_LINKOUT_URL = 'https://panelapp.genomicsengland.co.uk/panels/'
    GE_ZOOMA_DISEASE_MAPPING = '/tmp/zooma_disease_mapping.csv'
    GE_ZOOMA_DISEASE_MAPPING_NOT_HIGH_CONFIDENT = '/tmp/zooma_disease_mapping_low_confidence.csv'
    GE_PANEL_VERSION = 'v5.7'

    # IntoGEN
    INTOGEN_FILENAME = file_or_resource('intogen_opentargets.tsv')
    INTOGEN_EVIDENCE_FILENAME = HOME_DIR + '/otar001_intogen-18-12-2017.json'

    # mapping that we maintain in Zooma
    OMIM_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/xref_mappings/omim_to_efo.txt'
    ZOOMA_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/zooma/cttv_indications_3.txt'

    # mouse models
    MOUSEMODELS_PHENODIGM_SOLR = 'http://localhost:8983' # 'solrclouddev.sanger.ac.uk'
    # TODO remove refs to user directories
    MOUSEMODELS_CACHE_DIRECTORY = HOME_DIR + '/.phenodigmcache'
