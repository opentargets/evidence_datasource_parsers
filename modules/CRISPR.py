from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import requests
import python_jsonschema_objects as pjo

__copyright__ = "Copyright 2014-2019, Open Targets"
__credits__   = ["Michaela Spitzer"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "Open Targets Data Team"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

CRISPR_DATABASE_ID='Project Score'
CRISPR_VERSION='April 2019'

# A few genes do not have Ensembl IDs in the data file provided
CRISPR_SYMBOL_MAPPING={
    'CASC5': 'ENSG00000137812',
    'CIRH1A': 'ENSG00000141076',
    'EFTUD1': 'ENSG00000140598',
    'ENSG00000163660': 'ENSG00000163660',
    'KIAA0947': 'ENSG00000164151',
    'KIAA1432': 'ENSG00000107036',
    'NDNL2': 'ENSG00000185115',
    'SRPR': 'ENSG00000182934',
    'ZNF259': 'ENSG00000109917'
}

class CRISPR:
    def __init__(self, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):
        self.evidence_strings = list()
        self.symbols = {}
        self.logger = logging.getLogger(__name__)

        # Build JSON schema url from version
        self.schema_version = schema_version
        schema_url = "https://raw.githubusercontent.com/opentargets/json_schema/" + self.schema_version + "/draft4_schemas/opentargets.json"
        self.logger.info('Loading JSON schema at {}'.format(schema_url))

        # Initialize json builder based on the schema:
        try:
            r = requests.get(schema_url)
            r.raise_for_status()
            json_schema = r.json()
            self.builder = pjo.ObjectBuilder(json_schema)
            self.evidence_builder = self.builder.build_classes()
            self.schema_version = schema_version
        except requests.exceptions.HTTPError as e:
            self.logger.error('Invalid JSON schema version')
            raise e

    def process_crispr(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence()
        self.write_evidence()

    def build_evidence(self, filename1=Config.CRISPR_FILENAME1, filename2=Config.CRISPR_FILENAME2):

        now = datetime.datetime.now()

        # *** Build evidence.provenance_type object ***
        provenance_type = {
            'database': {
                'id': CRISPR_DATABASE_ID,
                'version': CRISPR_VERSION,
                'dbxref': {
                    'url': "https://score.depmap.sanger.ac.uk/",
                    'id': CRISPR_DATABASE_ID,
                    'version': CRISPR_VERSION,
                }
            }
        }

        # Build dictionary with publication info (pmid, method description, score type, min_score, max_score)
        CRISPR_METHOD_MAP = {}
        with open(filename2, 'r') as crispr_descriptions:

            for line in crispr_descriptions:
                # pmid	method	score_type
                (efo_id, crispr_method) = tuple(line.rstrip().split('\t'))
                CRISPR_METHOD_MAP[efo_id]={'method':crispr_method}

        with open(filename1, 'r') as crispr_input:
            n = 0
            for line in crispr_input:
                n +=1
                if n>1:
                    # pmid	gene_set_name	target_id	disease_id	disease_name	score
                    # 28892060	Intestine Key Driver Genes	CD53	EFO_0003767	Inflammatory bowel disease	0.083369979
                    # 28892060	Intestine Key Driver Genes	RHOH	EFO_0003767	Inflammatory bowel disease	0.11673251
                    # 28892060	Intestine Key Driver Genes	DOCK2	EFO_0003767	Inflammatory bowel disease	0.122212311
                    # 28892060	Intestine Key Driver Genes	FGR	EFO_0003767	Inflammatory bowel disease	0.13290268
                    (pmid, gene_set_name, target_name, disease_id, disease_name, score) = tuple(line.rstrip().split('\t'))


                    # *** Build evidence.resource_score object ***
                    resource_score = {
                        'type': "probability",
                        'method': {
                            'description': CRISPR_METHOD_MAP[disease_id]['method'],
                            'reference': "http://europepmc.org/abstract/MED/{0}".format(pmid)
                        },
                        'value': float(score)/100
                    }

                    # General properties of the evidence
                    validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                    access_level = "public"
                    type = "affected_pathway"
                    sourceID = "crispr"

                    if target_name in CRISPR_SYMBOL_MAPPING:
                        target_name=CRISPR_SYMBOL_MAPPING[target_name]

                    if target_name in self.symbols.values():

                        # *** Build target object ***
                        ensembl_gene_id = target_name
                        target_info = {
                            'id': "http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                            'target_name': target_name,
                            # TODO activity is a required field in target object, currently set as unknown
                            'activity': "http://identifiers.org/cttv.activity/unknown",
                            'target_type': "http://identifiers.org/cttv.target/gene_evidence"
                        }

                        # *** Build disease object ***
                        disease_info = {
                            'id': disease_id,
                            'name': disease_name
                        }

                        # *** Build evidence oject ***
                        evidence_info = {
                            'date_asserted': now.isoformat(),
                            'is_associated': True,
                            'evidence_codes': ["http://purl.obolibrary.org/obo/ECO_0000053"],
                            'provenance_type': provenance_type,
                            'resource_score': resource_score
                        }

                        # *** Build unique_association_field object ***
                        unique_association_fields = {
                            'pmid': "http://europepmc.org/abstract/MED/{0}".format(pmid),
                            ## TODO: Add gene_set_name to main part of evidence string
                            'gene_set': gene_set_name,
                            'target_id': target_info['id'],
                            'disease_id': disease_info['id']
                        }

                        try:
                            evidence_string = self.evidence_builder.Opentargets(
                                type = type,
                                access_level = access_level,
                                sourceID = sourceID,
                                evidence = evidence_info,
                                target = target_info,
                                disease = disease_info,
                                unique_association_fields = unique_association_fields,
                                validated_against_schema_version = validated_against_schema_version
                            )

                            # Append current evidence string
                            self.evidence_strings.append(evidence_string)
                        except:
                            self.logger.error('Evidence generation failed for the following row:\n{}'.format(line.rstrip()))
                            raise
                    else:
                        self.logger.error("%s is not found in Ensembl" % target_name)

            self.logger.error("%s evidence parsed"%(n-1))
            self.logger.error("%s evidence created"%len(self.evidence_strings))

    def write_evidence(self, filename=Config.CRISPR_EVIDENCE_FILENAME):
        self.logger.info(f"Writing CRISPR evidence strings to {Config.CRISPR_EVIDENCE_FILENAME}")

        with open(filename, 'w') as crispr_output:
            for evidence_string in self.evidence_strings:
                crispr_output.write(evidence_string.serialize() + "\n")

def main():
    CRISPR().process_crispr()

if __name__ == "__main__":
    main()

