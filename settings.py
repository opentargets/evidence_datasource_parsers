# General settings that all parsers can share

import os
import pkg_resources as res
from datetime import datetime

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
    # get filename and check if in getcwd then get from the package resources folder
    filename = os.path.expanduser(fname)

    resource_package = __name__
    resource_path = '/'.join(('resources', filename))

    if filename is not None:
        abs_filename = os.path.join(os.path.abspath(os.getcwd()), filename) \
                       if not os.path.isabs(filename) else filename

        return abs_filename if os.path.isfile(abs_filename) \
            else res.resource_filename(resource_package, resource_path)

class Config:
    # shared settings

    # schema version
    VALIDATED_AGAINST_SCHEMA_VERSION = '1.2.8'

    # Ontologies
    EFO_URL = 'https://github.com/EBISPOT/efo/raw/v2018-01-15/efo.obo'
    HP_URL = 'http://purl.obolibrary.org/obo/hp.obo'

    # HGNC
    GENES_HGNC = 'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'

    # UKBIOBANK
    UKBIOBANK_FILENAME = file_or_resource('ukbiobank.txt')
    UKBIOBANK_EVIDENCE_FILENAME = 'ukbiobank-30-04-2018.json'

    # SYSBIO
    SYSBIO_FILENAME1 = file_or_resource('sysbio_evidence-31-01-2019.tsv')
    SYSBIO_FILENAME2 = file_or_resource('sysbio_publication_info_nov2018.tsv')
    SYSBIO_EVIDENCE_FILENAME = 'sysbio-29-01-2019.json'

    # mapping that we maintain in Zooma
    OMIM_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/xref_mappings/omim_to_efo.txt'
    ZOOMA_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/zooma/cttv_indications_3.txt'

    # Configuration for genetics portal evidences:
    ACTIVITY_URL = 'http://identifiers.org/cttv.activity'
    TARGET_URL = 'http://identifiers.org/ensembl'
    TARGET_TYPE_URL = 'http://identifiers.org/cttv.target'
    LITERATURE_URL = 'http://europepmc.org/abstract/MED'
    GENETICS_PORTAL_URL = 'https://genetics.opentargets.org'
    DISEASE_URL = 'http://www.ebi.ac.uk/efo'
    CONSEQUENCE_URL = 'http://purl.obolibrary.org/obo'

    # Evidence codes:
    EVIDENCE_CODE_INFERENCE = 'http://purl.obolibrary.org/obo/ECO_0000362' # computational inference
    EVIDENCE_CODE_EVIDENCE_TYPE = 'http://identifiers.org/eco/GWAS' # GWAS data type.
    EVIDENCE_CODE_SOURCE = 'http://identifiers.org/eco/locus_to_gene_pipeline' # variant to gene derived from l2g pipeline
