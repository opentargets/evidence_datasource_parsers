import os
import pkg_resources as res
import configparser

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
    PHEWAS_CATALOG_LOCN = 'https://storage.googleapis.com/phewas-catalog'
    GENES_HGNC =  'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'
    PHEWAS_CATALOG_JSON = '../phewas.json'

    INTOGEN_FILENAME = file_or_resource('intogen_opentargets.tsv')
    INTOGEN_EVIDENCE_FILENAME = '/tmp/intogen.json'

    ONTOLOGY_CONFIG = configparser.ConfigParser()
    ONTOLOGY_CONFIG.read(file_or_resource('ontology_config.ini'))

    # mapping that we maintain in Zooma
    OMIM_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/xref_mappings/omim_to_efo.txt'
    ZOOMA_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/zooma/cttv_indications_3.txt'

    # GE Pipeline

    GE_EVIDENCE_STRING = '/tmp/genomics_england_evidence_string.json'
    GE_LINKOUT_URL = 'https://panelapp.genomicsengland.co.uk/panels/'
    GE_ZOOMA_DISEASE_MAPPING = '/tmp/zooma_disease_mapping.csv'
    GE_ZOOMA_DISEASE_MAPPING_NOT_HIGH_CONFIDENT = '/tmp/zooma_disease_mapping_low_confidence.csv'

    VALIDATED_AGAINST_SCHEMA_VERSION = '1.2.7'

    MOUSEMODELS_CACHE_DIRECTORY = '/Users/otvisitor/.phenodigmcache'