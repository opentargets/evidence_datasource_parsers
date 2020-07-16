from settings import Config
from common.HGNCParser import GeneParser
from common.RareDiseasesUtils import RareDiseaseMapper

import python_jsonschema_objects as pjo

import logging
import pandas as pd
import gzip
import requests
import datetime
import argparse

__copyright__  = "Copyright 2014-2020, Open Targets"
__credits__    = ["Asier Gonzalez" ]
__license__    = "Apache 2.0"
__version__    = "0.0.1"
__maintainer__ = "Open Targets Data Team"
__email__      = ["data@opentargets.org"]
__status__     = "Prototype"

ClinGen_classification2score = {
    "Definitive": 1,
    "Strong": 0.75,
    "Moderate": 0.5,
    "Limited": 0.25,
    "Disputed": 0.1,
    "Refuted": 0.05,
    "No Reported Evidence": 0.01,
}

class ClinGen(RareDiseaseMapper):
    def __init__(self, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):
        super(ClinGen, self).__init__()
        self.genes = None
        self.evidence_strings = list()

        self._logger = logging.getLogger(__name__)

        # Build JSON schema url from version
        self.schema_version = schema_version
        schema_url = "https://raw.githubusercontent.com/opentargets/json_schema/" + self.schema_version + "/draft4_schemas/opentargets.json"
        self._logger.info('Loading JSON schema at {}'.format(schema_url))

        # Initialize json builder based on the schema:
        try:
            r = requests.get(schema_url)
            r.raise_for_status()
            json_schema = r.json()
            self.builder = pjo.ObjectBuilder(json_schema)
            self.evidence_builder = self.builder.build_classes()
            self.schema_version = schema_version
        except requests.exceptions.HTTPError as e:
            self._logger.error('Invalid JSON schema version')
            raise e


    def process_gene_validity_curations(self, in_filename, out_filename):

        self.get_omim_to_efo_mappings()

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        self.generate_evidence_strings(in_filename)

        # Save results to file
        self.write_evidence_strings(out_filename)

    def generate_evidence_strings(self, filename):

        total_efo = 0

        # When reading csv file skip header lines that don't contain column names
        # CLINGEN GENE VALIDITY CURATIONS
        # FILE CREATED: 2020-07-16
        # WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
        # +++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
        # GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE
        # +++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++

        gene_validity_curation_df = pd.read_csv(filename, skiprows= [0,1,2,3,5], quotechar='"')
        for index, row in gene_validity_curation_df.iterrows():
            print("{} - {}".format(row["GENE SYMBOL"], row["DISEASE LABEL"]))

            gene_symbol = row["GENE SYMBOL"]
            disease_name = row["DISEASE LABEL"]
            disease_id = row["DISEASE ID"]
            mode_of_inheritancr = row["MOI"]
            classification = row["CLASSIFICATION"]
            date = row["CLASSIFICATION DATE"]

            gene_symbol.rstrip()

            if gene_symbol in self.genes:
                # Map gene symbol to ensembl
                target = self.genes[gene_symbol]
                ensembl_iri = "http://identifiers.org/ensembl/" + target


                type = "genetic_literature"

                provenance_type = {
                    'database' : {
                        'id' : "ClinGen - Gene Validity Curations",
                        'version' : '2020.07.16',
                        'dbxref' : {
                            'url': "https://search.clinicalgenome.org/kb/gene-validity",
                            'id' : "ClinGen - Gene Validity Curations",
                            'version' : "2020.07.16"

                        }
                    }
                }

                # *** General properties ***
                access_level = "public"
                sourceID = "clingen"
                validated_against_schema_version = self.schema_version

                # *** Target info ***
                target = {
                    'id' : ensembl_iri,
                    'activity' : "http://identifiers.org/cttv.activity/unknown",
                    'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                    'target_name' : gene_symbol
                }
                # http://www.ontobee.org/ontology/ECO?iri=http://purl.obolibrary.org/obo/ECO_0000204 -- An evidence type that is based on an assertion by the author of a paper, which is read by a curator.

                # *** Disease info ***
                disease_info = {
                    'id' : "http://purl.obolibrary.org/obo/MONDO_0019181" + disease_id,
                    #'name' : disease['efo_label'],
                    'source_name' : disease_name
                }
                # *** Evidence info ***
                # Score based on disease confidence/ classification
                if classification in ClinGen_classification2score:
                    score = ClinGen_classification2score[classification]
                else:
                    self.logger.error('{} is not a recognised ClinGen classification, assigning an score of 0'.format(classification))
                    score = 0
                resource_score = {
                    'type': "probability",
                    'value': score
                }

                # Linkout
                linkout = [
                    {
                        'url' : 'http://www.ebi.ac.uk/gene2phenotype/search?panel=ALL&search_term=%s' % (gene_symbol,),
                        'nice_name' : 'Gene2Phenotype%s' % (gene_symbol)
                    }
                ]

                evidence = {
                    'is_associated' : True,
                    'confidence' : confidence,
                    'allelic_requirement' : allelic_requirement,
                    'mutation_consequence' : mutation_consequence,
                    'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000204"],
                    'provenance_type' : provenance_type,
                    'date_asserted' : date,
                    'resource_score' : resource_score,
                    'urls' : linkout
                }
                # *** unique_association_fields ***
                unique_association_fields = {
                    'target_id' : ensembl_iri,
                    'original_disease_label' : disease_name,
                    'disease_id' : disease['efo_uri'],
                    'gene_panel': panel
                }


                try:
                    evidence = self.evidence_builder.Opentargets(
                        type = type,
                        access_level = access_level,
                        sourceID = sourceID,
                        evidence = evidence,
                        target = target,
                        disease = disease_info,
                        unique_association_fields = unique_association_fields,
                        validated_against_schema_version = validated_against_schema_version
                    )
                    self.evidence_strings.append(evidence)
                except:
                    self._logger.warning('Evidence generation failed for row: {}'.format(c))
                    raise

    def write_evidence_strings(self, filename):
        self._logger.info("Writing Gene2Phenotype evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self._logger.info(evidence_string['disease']['id'])
                tp_file.write(evidence_string.serialize() + "\n")
        tp_file.close()


def main():

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-i', '--input_file',
                        help='Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file',
                        type=str, required=True)
    parser.add_argument('-s', '--schema_version',
                        help='JSON schema version to use, e.g. 1.6.8. It must be branch or a tag available in https://github.com/opentargets/json_schema',
                        type=str, required=True)

    args = parser.parse_args()
    # Get parameters
    infile = args.input_file
    outfile = args.output_file
    schema_version = args.schema_version

    clingen = ClinGen(schema_version=schema_version)
    clingen.process_gene_validity_curations(infile, outfile)

if __name__ == "__main__":
    main()