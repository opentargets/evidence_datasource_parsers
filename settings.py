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

    GOOGLE_DEFAULT_PROJECT = 'open-targets'
    GOOGLE_BUCKET_EVIDENCE_INPUT = 'otar000-evidence_input'

    # Ontologies
    EFO_URL = 'https://github.com/EBISPOT/efo/raw/v2018-01-15/efo.obo'
    HP_URL = 'http://purl.obolibrary.org/obo/hp.obo'

    # HGNC
    GENES_HGNC = 'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'

    # PROGENY
    PROGENY_FILENAME = file_or_resource('progeny_normalVStumor_opentargets.txt')
    PROGENY_EVIDENCE_FILENAME = 'progeny-20-05-2018.json'

    # UKBIOBANK
    UKBIOBANK_FILENAME = file_or_resource('ukbiobank.txt')
    UKBIOBANK_EVIDENCE_FILENAME = 'ukbiobank-30-04-2018.json'

    # SYSBIO
    SYSBIO_FILENAME1 = file_or_resource('sysbio_evidence-31-01-2019.tsv')
    SYSBIO_FILENAME2 = file_or_resource('sysbio_publication_info_nov2018.tsv')
    SYSBIO_EVIDENCE_FILENAME = 'sysbio-29-01-2019.json'

    # CRISPR
    CRISPR_FILENAME1 = file_or_resource('crispr_evidence.tsv')
    CRISPR_FILENAME2 = file_or_resource('crispr_descriptions.tsv')
    CRISPR_EVIDENCE_FILENAME = 'crispr-21-08-2019.json'

    # Gene2Phenotype
    #G2P_FILENAME = 'DDG2P.csv.gz'
    G2P_DD_FILENAME = file_or_resource('DDG2P_2_4_2020.csv.gz')
    G2P_eye_FILENAME = file_or_resource('EyeG2P_26_3_2020.csv.gz')
    G2P_skin_FILENAME = file_or_resource('SkinG2P_26_3_2020.csv.gz')
    G2P_cancer_FILENAME = file_or_resource('CancerG2P_26_3_2020.csv.gz')
    G2P_EVIDENCE_FILENAME = 'gene2phenotype-19-08-2019.json'

    # IntoGEN
    INTOGEN_DRIVER_GENES_FILENAME = file_or_resource('intogen_Compendium_Cancer_Genes.tsv')
    INTOGEN_EVIDENCE_FILENAME = 'intogen-02-02-2020.json'
    INTOGEN_CANCER2EFO_MAPPING_FILENAME = file_or_resource('intogen_cancer2EFO_mapping.tsv')
    INTOGEN_COHORTS = file_or_resource('intogen_cohorts.tsv')

    # mapping that we maintain in Zooma
    OMIM_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/xref_mappings/omim_to_efo.txt'
    ZOOMA_TO_EFO_MAP_URL = 'https://raw.githubusercontent.com/opentargets/platform_semantic/master/resources/zooma/cttv_indications_3.txt'

    # mouse models (Phenodigm)
    # used to be 'http://localhost:8983' # 'solrclouddev.sanger.ac.uk'
    MOUSEMODELS_PHENODIGM_SOLR = 'http://www.ebi.ac.uk/mi/impc'
    # write to the cloud direcly
    MOUSEMODELS_CACHE_DIRECTORY = 'PhenoDigm/phenodigmcache'
    MOUSEMODELS_EVIDENCE_FILENAME = f"phenodigm-{datetime.today().strftime('%Y-%m-%d')}.json.gz"

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
