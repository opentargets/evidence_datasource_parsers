import argparse
import json
import logging

import pandas as pd
import ontoma

class ClinGen():

    def __init__(self):

        self.evidence_strings = list()
        self.unmapped_diseases = set()

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
        gene_validity_curation_df = pd.read_csv(filename, skiprows=[0, 1, 2, 3, 5], quotechar='"')

        for index, row in gene_validity_curation_df.iterrows():
            logging.info("{} - {}".format(row["GENE SYMBOL"], row["DISEASE LABEL"]))

            gene_symbol = row["GENE SYMBOL"]
            disease_name = row["DISEASE LABEL"]
            disease_id = row["DISEASE ID (MONDO)"]
            mode_of_inheritance = row["MOI"]
            classification = row["CLASSIFICATION"]
            expert_panel_name = row["GCEP"]
            report_url = row["ONLINE REPORT"]

            gene_symbol.rstrip()

            # Default value for efo mappings: 
            efo_mappings = [ None ]

            # Check that disease id exists in EFO or find equivalent term
            # If id is not found in EFO OBO file the function returns None
            if self.ontoma.get_efo_label(disease_id):
                disease_label = self.ontoma.get_efo_label(disease_id)
                # Create list of single disease to mimic what is returned by next step
                efo_mappings = [{'id': disease_id, 'name': disease_label}]
            
            elif self.ontoma.get_efo_from_xref(disease_id):
                efo_mappings = self.ontoma.get_efo_from_xref(disease_id)
                logging.info("{} mapped to {} EFO ids based on xrefs.".format(disease_id, len(efo_mappings)))
            
            else:
                # Search disease label using OnToma and accept perfect matches
                ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
                if ontoma_mapping:
                    if ontoma_mapping['action'] is None:
                        efo_mappings = [{'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}]
                    else:
                        # OnToma fuzzy match ignored
                        logging.warning("Fuzzy match from OnToma ignored. Request EFO team to import {} - {}".format(disease_name, disease_id))
                        # Record the unmapped disease
                        self.unmapped_diseases.add((disease_id, disease_name))
                        
                else:
                    # MONDO id could not be found in EFO. Log it and continue
                    logging.warning("{} - {} could not be mapped to any EFO id. Skipping it, it should be checked with the EFO team".format(disease_name, disease_id))
                    
                    # Record the unmapped disease
                    self.unmapped_diseases.add((disease_id, disease_name))

            for efo_mapping in efo_mappings:

                evidence = {
                    'datasourceId': 'clingen',
                    'datatypeId': 'genetic_literature',
                    'targetFromSourceId': gene_symbol,
                    'diseaseFromSource': disease_name,
                    'diseaseFromSourceId': disease_id,
                    'allelicRequirements': [mode_of_inheritance],
                    'confidence': classification,
                    'studyId': expert_panel_name,
                    'urls': {'url': report_url}
                }

                # diseaseFromSourceMappedId might be None - in that case the key is skipped:
                if efo_mapping:
                    evidence['diseaseFromSourceMappedId'] = ontoma.interface.make_uri(efo_mapping['id']).split("/")[-1]

                self.evidence_strings.append(evidence)

    def write_evidence_strings(self, filename):
        logging.info("Writing ClinGen evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            for evidence_string in self.evidence_strings:
                json.dump(evidence_string, tp_file)
                tp_file.write("\n")

    def write_unmapped_diseases(self, filename):
        logging.info("Writing ClinGen diseases not mapped to EFO to %s", filename)
        with open(filename, 'w') as unmapped_diseases_file:
            unmapped_diseases_file.write("disease_id\tdisease_name\n")
            for unmapped_disease in self.unmapped_diseases:
                unmapped_diseases_file.write(unmapped_disease[0] + "\t" + unmapped_disease[1] + "\n")


def main(infile, outfile, unmapped_diseases_file):

    clingen = ClinGen()
    clingen.process_gene_validity_curations(infile, outfile, unmapped_diseases_file)


if __name__ == "__main__":

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
    parser.add_argument('-l', '--log_file', type=str,
                        help='Optional filename to redirect the logs into.')

    args = parser.parse_args()

    # Get parameters
    infile = args.input_file
    outfile = args.output_file
    unmapped_diseases_file = args.unmapped_diseases_file
    log_file = args.log_file

    # Initialize the logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    main(infile, outfile, unmapped_diseases_file)
