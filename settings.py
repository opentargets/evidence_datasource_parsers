import os
import pkg_resources as res
import configparser
from pathlib import Path
from envparse import env, ConfigurationError

def read_option(option, cast=None,
                **kwargs):

    try:
        default_value = kwargs.pop('default')
    except KeyError:
        default_value = None

    try:
        # reading the environment variable with envparse
        return env(option, cast=cast, **kwargs)
    except ConfigurationError:
       return default_value


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

    HOME_DIR = str(Path.home())

    # schema version
    VALIDATED_AGAINST_SCHEMA_VERSION = '1.2.7'

    GOOGLE_DEFAULT_PROJECT = 'open-targets'
    GOOGLE_BUCKET_EVIDENCE_INPUT = 'otar000-evidence_input'

    #Ontologies
    EFO_URL = 'https://github.com/EBISPOT/efo/raw/v2018-01-15/efo.obo'
    HP_URL = 'http://purl.obolibrary.org/obo/hp.obo'

    # HGNC
    GENES_HGNC =  'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'

    # SLAPEnrich
    SLAPENRICH_FILENAME = file_or_resource('slapenrich_opentargets.tsv')
    SLAPENRICH_EVIDENCE_FILENAME = HOME_DIR + '/otar001_slapenrich-19-12-2017.json'

    # Gene2Phenotype
    #G2P_FILENAME = file_or_resource('DDG2P_14_5_2017.csv.gz')
    G2P_FILENAME = 'DDG2P_14_5_2017.csv.gz'
    G2P_EVIDENCE_FILENAME = 'gene2phenotype.json'

    # Genomics England
    GE_EVIDENCE_STRING = HOME_DIR + '/otar001_genomics_england-18-12-2017.json'
    GE_LINKOUT_URL = 'https://panelapp.genomicsengland.co.uk/panels/'
    GE_ZOOMA_DISEASE_MAPPING = '/tmp/zooma_disease_mapping.csv'
    GE_ZOOMA_DISEASE_MAPPING_NOT_HIGH_CONFIDENT = '/tmp/zooma_disease_mapping_low_confidence.csv'

    # IntoGEN
    INTOGEN_FILENAME = file_or_resource('intogen_opentargets.tsv')
    INTOGEN_EVIDENCE_FILENAME = HOME_DIR + '/otar001_intogen-18-12-2017.json'

    # Phewas
    # PHEWAS_CATALOG_URL = 'https://phewascatalog.org/files/phewas-catalog.csv.zip'
    PHEWAS_CATALOG_URL = 'https://storage.googleapis.com/otar000-evidence_input/PheWAScatalog/phewas-catalog.csv'
    PHEWAS_PHECODE_MAP_URL = 'https://phewascatalog.org/files/phecode_icd9_map_unrolled.csv.zip'
    PHEWAS_CATALOG_LOCN = 'https://storage.googleapis.com/phewas-catalog'
    PHEWAS_CATALOG_JSON = HOME_DIR + '/otar001_phewas_catalog-19-12-2017.json'

    ONTOLOGY_CONFIG = configparser.ConfigParser()
    ONTOLOGY_CONFIG.read(file_or_resource('ontology_config.ini'))

    # mapping that we maintain in Zooma
    OMIM_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/xref_mappings/omim_to_efo.txt'
    ZOOMA_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/zooma/cttv_indications_3.txt'

    # mouse models
    MOUSEMODELS_PHENODIGM_SOLR = 'http://localhost:8983' # 'solrclouddev.sanger.ac.uk'
    # TODO remove refs to user directories
    MOUSEMODELS_CACHE_DIRECTORY =  HOME_DIR + '/.phenodigmcache'

    # MONGO_URL = read_option('MONGO_URL', cast=str, default='')
    # MONGO_DB = read_option('MONGO_DB', cast=str, default='')
    # MONGO_TABLE = read_option('MONGO_TABLE', cast=str, default='')
    # MONGO_USER = read_option('MONGO_USER', cast=str, default='')
    # MONGO_PWD = read_option('MONGO_PWD', cast=str, default='')
