from settings import Config
from common.HGNCParser import GeneParser
import ontoma

import python_jsonschema_objects as pjo

import logging
import pandas as pd
import requests
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
    "Strong": 1,
    "Moderate": 0.5,
    "Limited": 0.01,
    "Disputed": 0.01,
    "Refuted": 0.01,
    "No Reported Evidence": 0.01,
}

class ClinGen():
    def __init__(self, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):
        self.genes = None
        self.evidence_strings = list()
        self.unmapped_diseases = set()

        # Configure logging
        # Create logger
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(logging.INFO)

        # Create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)

        # Create formatter
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

        # Add formatter to ch
        ch.setFormatter(formatter)

        # Add ch to handler
        self._logger.addHandler(ch)

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

        # Create OnToma object
        self.ontoma = ontoma.interface.OnToma()


    def process_gene_validity_curations(self, in_filename, out_filename, unmapped_diseases_filename):

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        self.generate_evidence_strings(in_filename)

        # Save results to file
        self.write_evidence_strings(out_filename)

        # If selected, write unmapped diseases to file
        if unmapped_diseases_filename:
            self.write_unmapped_diseases(unmapped_diseases_filename)

    def generate_evidence_strings(self, filename):

        # When reading csv file first extract date from second row and the skip header lines that don't contain column names
        # CLINGEN GENE VALIDITY CURATIONS
        # FILE CREATED: 2020-07-16
        # WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
        # +++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
        # GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE
        # +++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++

        # Read first two rows of file to extract date
        with open(filename) as f:
            f.readline() # Don't do anything with first row
            second_row = f.readline()
            file_created_date = second_row.split(": ")[1].strip() # Extract date and clean it


        gene_validity_curation_df = pd.read_csv(filename, skiprows= [0,1,2,3,5], quotechar='"')
        for index, row in gene_validity_curation_df.iterrows():
            print("{} - {}".format(row["GENE SYMBOL"], row["DISEASE LABEL"]))

            gene_symbol = row["GENE SYMBOL"]
            disease_name = row["DISEASE LABEL"]
            disease_id = row["DISEASE ID (MONDO)"]
            mode_of_inheritancr = row["MOI"]
            classification = row["CLASSIFICATION"]
            date = row["CLASSIFICATION DATE"]
            report_url = row["ONLINE REPORT"]

            gene_symbol.rstrip()

            if gene_symbol in self.genes:
                # Map gene symbol to ensembl
                target = self.genes[gene_symbol]
                ensembl_iri = "http://identifiers.org/ensembl/" + target

                # Check that disease id exists in EFO or find equivalent term
                # If id is not found in EFO OBO file the function returns None
                if self.ontoma.get_efo_label(disease_id):
                    disease_label = self.ontoma.get_efo_label(disease_id)
                    # Create list of single disease to mimic what is returned by next step
                    efo_mappings = [{'id': disease_id, 'name': disease_label}]
                elif self.ontoma.get_efo_from_xref(disease_id):
                    efo_mappings = self.ontoma.get_efo_from_xref(disease_id)
                    self._logger.info("{} mapped to {} EFO ids based on xrefs.".format(disease_id, len(efo_mappings)))
                else:
                    # Search disease label using OnToma and accept perfect matches
                    ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
                    if ontoma_mapping:
                        if ontoma_mapping['action'] is None:
                            efo_mappings = [{'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}]
                        else:
                            # OnToma fuzzy match ignored
                            self._logger.info("Fuzzy match from OnToma ignored. Request EFO team to import {} - {}".format(disease_name, disease_id))
                            # Record the unmapped disease
                            self.unmapped_diseases.add((disease_id, disease_name))
                            continue
                    else:
                        # MONDO id could not be found in EFO. Log it and continue
                        self._logger.info("{} - {} could not be mapped to any EFO id. Skipping it, it should be checked with the EFO team".format(disease_name, disease_id))
                        # Record the unmapped disease
                        self.unmapped_diseases.add((disease_id, disease_name))
                        continue

                for efo_mapping in efo_mappings:

                    # *** Disease info ***
                    disease_info = {
                        'id': ontoma.interface.make_uri(efo_mapping['id']),
                        'name' : efo_mapping['name'],
                        'source_name': disease_name
                    }

                    type = "genetic_literature"

                    provenance_type = {
                        'database' : {
                            'id' : "ClinGen - Gene Validity Curations",
                            'version' : file_created_date,
                            'dbxref' : {
                                'url': "https://search.clinicalgenome.org/kb/gene-validity",
                                'id' : "ClinGen - Gene Validity Curations",
                                'version' : file_created_date

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
                            'url' : report_url,
                            'nice_name' : 'Gene Validity Curations: {} - {} report'.format(gene_symbol, disease_name)
                        }
                    ]

                    evidence = {
                        'is_associated' : True,
                        'confidence' : classification,
                        'allelic_requirement' : mode_of_inheritancr,
                        'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000204"],
                        'provenance_type' : provenance_type,
                        'date_asserted' : date,
                        'resource_score' : resource_score,
                        'urls' : linkout
                    }
                    # *** unique_association_fields ***
                    unique_association_fields = {
                        'target_id' : ensembl_iri,
                        'disease_id' : disease_id,
                        'report_url' : report_url
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
                        self._logger.warning('Evidence generation failed for row: {}'.format(index))
                        raise

    def write_evidence_strings(self, filename):
        self._logger.info("Writing ClinGen evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            for evidence_string in self.evidence_strings:
                tp_file.write(evidence_string.serialize() + "\n")

    def write_unmapped_diseases(self, filename):
        self._logger.info("Writing ClinGen diseases not mapped to EFO to %s", filename)
        with open(filename, 'w') as unmapped_diseases_file:
            unmapped_diseases_file.write("disease_id\tdisease_name\n")
            for unmapped_disease in self.unmapped_diseases:
                unmapped_diseases_file.write(unmapped_disease[0] + "\t" + unmapped_disease[1] + "\n")


def main():

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse ClinGen gene-disease associations from Gene Validity Curations')
    parser.add_argument('-i', '--input_file',
                        help='Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file',
                        type=str, required=True)
    parser.add_argument('-s', '--schema_version',
                        help='JSON schema version to use, e.g. 1.6.8. It must be branch or a tag available in https://github.com/opentargets/json_schema',
                        type=str, required=True)
    parser.add_argument('-u', '--unmapped_diseases_file',
                        help='If specified, the diseases not mapped to EFO will be stored in this file',
                        type=str, default=False)

    args = parser.parse_args()
    # Get parameters
    infile = args.input_file
    outfile = args.output_file
    schema_version = args.schema_version
    unmapped_diseases_file = args.unmapped_diseases_file

    clingen = ClinGen(schema_version=schema_version)
    clingen.process_gene_validity_curations(infile, outfile, unmapped_diseases_file)

if __name__ == "__main__":
    main()