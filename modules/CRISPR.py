from settings import Config
from common.HGNCParser import GeneParser

import logging
import datetime
import requests
import pandas as pd
import python_jsonschema_objects as pjo
import argparse

__copyright__ = "Copyright 2014-2020, Open Targets"
__credits__   = ["Michaela Spitzer", "Asier Gonzalez"]
__license__   = "Apache 2.0"
__version__   = "1.3.0"
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

    def process_crispr(self, in_filename, desc_filename, cell_filename, out_filename):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence(in_filename, desc_filename, cell_filename)
        self.write_evidence(out_filename)

    def build_evidence(self, input_filename=Config.CRISPR_FILENAME1, description_filename=Config.CRISPR_FILENAME2, cell_line_filename=None):

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

        # Get cell lines per cancer type
        crispr_desc_df = pd.read_csv(description_filename ,sep='\t')
        crispr_cells_df = pd.read_csv(cell_line_filename, sep='\t')
        # Join using cancer type
        joined_cancer_type_df = crispr_desc_df.join(crispr_cells_df.set_index('Cancer Type'),
                                                    on='tissue_or_cancer_type', how='inner').drop('Tissue', axis=1)
        # Join using tissue
        joined_tissue_df = crispr_desc_df.join(crispr_cells_df.set_index('Tissue'), on='tissue_or_cancer_type',
                                               how='inner').drop('Cancer Type', axis=1)
        # Concatenate two data frames
        concat_df = pd.concat([joined_cancer_type_df, joined_tissue_df ], ignore_index=True)

        # Aggregate cell lines in Name column
        joined_df = concat_df.groupby(['efo_id', 'tissue_or_cancer_type', 'method']).agg(
            {'Name': lambda x: list(x)}).reset_index()
        joined_df = joined_df.set_index('efo_id')

        with open(input_filename, 'r') as crispr_input:
            n = 0
            for line in crispr_input:
                n +=1
                if n>1:
                    # pmid	gene_set_name	target_id	disease_id	disease_name	score
                    # 30971826	Project Score: Prioritisation of oncology therapeutic targets using CRISPR-Cas9 screening	ENSG00000110092	http://www.ebi.ac.uk/efo/EFO_0000174	Ewing sarcoma	79.7813
                    # 30971826	Project Score: Prioritisation of oncology therapeutic targets using CRISPR-Cas9 screening	ENSG00000130725	http://www.ebi.ac.uk/efo/EFO_0000174	Ewing sarcoma	58.5
                    # 30971826	Project Score: Prioritisation of oncology therapeutic targets using CRISPR-Cas9 screening	ENSG00000111142	http://www.ebi.ac.uk/efo/EFO_0000174	Ewing sarcoma	57.9375
                    # 30971826	Project Score: Prioritisation of oncology therapeutic targets using CRISPR-Cas9 screening	ENSG00000152234	http://www.ebi.ac.uk/efo/EFO_0000174	Ewing sarcoma	54.125
                    # 30971826	Project Score: Prioritisation of oncology therapeutic targets using CRISPR-Cas9 screening	ENSG00000165501	http://www.ebi.ac.uk/efo/EFO_0000174	Ewing sarcoma	52.5
                    (pmid, gene_set_name, target_name, disease_id, disease_name, score) = tuple(line.rstrip().split('\t'))


                    # *** Build evidence.resource_score object ***
                    resource_score = {
                        'type': "probability",
                        'method': {
                            'description': joined_df.loc[disease_id]['method'],
                            'reference': "http://europepmc.org/abstract/MED/{0}".format(pmid)
                        },
                        'value': float(score)/100
                    }

                    # General properties of the evidence
                    validated_against_schema_version = self.schema_version
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
                            'target_type': "http://identifiers.org/cttv.target/gene_evidence"
                        }

                        # *** Build disease object ***
                        disease_info = {
                            'id': disease_id,
                            'name': disease_name,
                            'source_name': joined_df.loc[disease_id]['tissue_or_cancer_type'] #NOTE: It may be a tissue or a cancer type
                        }

                        # *** Build evidence oject ***
                        evidence_info = {
                            'date_asserted': now.isoformat(),
                            'is_associated': True,
                            'evidence_codes': ["http://purl.obolibrary.org/obo/ECO_0000053"],
                            'provenance_type': provenance_type,
                            'resource_score': resource_score,
                            'experiment_overview': gene_set_name, #TODO: Add this to the schema
                            'cell_lines': joined_df.loc[disease_id]['Name'] #TODO: Add this to the schema
                        }

                        # *** Build unique_association_field object ***
                        unique_association_fields = {
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

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse essential cancer genes identified using CRISPR assays and prioritised in Project Score')
    parser.add_argument('-d', '--descriptions_file',
                        help='Name of tsv file with the description of the method per cancer type',
                        type=str, default=Config.CRISPR_FILENAME2)
    parser.add_argument('-i', '--input_file',
                        help='Name of tsv file with the priority score',
                        type=str, default=Config.CRISPR_FILENAME1)
    parser.add_argument('-c', '--cell_types_file',
                        help='Name of tsv file with cell line names per cancer type',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file',
                        type=str, default=Config.CRISPR_EVIDENCE_FILENAME)
    parser.add_argument('-s', '--schema_version',
                        help='JSON schema version to use, e.g. 1.7.5. It must be branch or a tag available in https://github.com/opentargets/json_schema',
                        type=str, required=True)

    args = parser.parse_args()
    # Get parameters
    desc_file = args.descriptions_file
    infile = args.input_file
    cell_file = args.cell_types_file
    outfile = args.output_file
    schema_version = args.schema_version

    crispr_parser = CRISPR(schema_version)
    crispr_parser.process_crispr(infile, desc_file, cell_file, outfile)

if __name__ == "__main__":
    main()

