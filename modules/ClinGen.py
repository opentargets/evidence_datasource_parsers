import ontoma

import logging
import pandas as pd
import argparse
import json


class ClinGen():
    def __init__(self):

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

        # Create OnToma object
        self.ontoma = ontoma.interface.OnToma()


    def process_gene_validity_curations(self, in_filename, out_filename, unmapped_diseases_filename):

        self.generate_evidence_strings(in_filename)

        # Save results to file
        self.write_evidence_strings(out_filename)

        # If selected, write unmapped diseases to file
        if unmapped_diseases_filename:
            self.write_unmapped_diseases(unmapped_diseases_filename)

    def generate_evidence_strings(self, filename):

        # When reading csv file skip header lines that don't contain column names
        gene_validity_curation_df = pd.read_csv(filename, skiprows= [0,1,2,3,5], quotechar='"')
        for index, row in gene_validity_curation_df.iterrows():
            print("{} - {}".format(row["GENE SYMBOL"], row["DISEASE LABEL"]))

            gene_symbol = row["GENE SYMBOL"]
            disease_name = row["DISEASE LABEL"]
            disease_id = row["DISEASE ID (MONDO)"]
            mode_of_inheritance = row["MOI"]
            classification = row["CLASSIFICATION"]
            expert_panel_name = row["GCEP"]
            report_url = row["ONLINE REPORT"]

            gene_symbol.rstrip()

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

                linkout = [
                    {
                        'url' : report_url
                    }
                ]

                evidence = {
                    'datasourceId' : 'clingen',
                    'datatypeId' : 'genetic_literature',
                    'targetFromSourceId' : gene_symbol,
                    'diseaseFromSource': disease_name,
                    'diseaseFromSourceId': disease_id,
                    'diseaseFromSourceMappedId': ontoma.interface.make_uri(efo_mapping['id']).split("/")[-1],
                    'allelicRequirements': [mode_of_inheritance],
                    'confidence': classification,
                    'studyId': expert_panel_name,
                    'urls': linkout
                }

                self.evidence_strings.append(evidence)

    def write_evidence_strings(self, filename):
        self._logger.info("Writing ClinGen evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            for evidence_string in self.evidence_strings:
                json.dump(evidence_string, tp_file)
                tp_file.write("\n")

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
    parser.add_argument('-u', '--unmapped_diseases_file',
                        help='If specified, the diseases not mapped to EFO will be stored in this file',
                        type=str, default=False)

    args = parser.parse_args()
    # Get parameters
    infile = args.input_file
    outfile = args.output_file
    unmapped_diseases_file = args.unmapped_diseases_file

    clingen = ClinGen()
    clingen.process_gene_validity_curations(infile, outfile, unmapped_diseases_file)

if __name__ == "__main__":
    main()