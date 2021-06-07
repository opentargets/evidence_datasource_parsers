import logging
import csv
import gzip
import argparse
import json
import sys

import ontoma


G2P_mutationCsq2functionalCsq = {
    'loss of function': 'SO_0002054',  # loss_of_function_variant
    'all missense/in frame': 'SO_0001650',  # inframe_variant
    'uncertain': 'SO_0002220',  # function_uncertain_variant
    'activating': 'SO_0002053',  # gain_of_function_variant
    'dominant negative': 'SO_0002052',  # dominant_negative_variant
    '': None,
    'gain of function': 'SO_0002053',  # gain_of_function_variant
    'cis-regulatory or promotor mutation': 'SO_0001566',  # regulatory_region_variant
    '5_prime or 3_prime UTR mutation': 'SO_0001622',  # UTR_variant
    'increased gene dosage': 'SO_0001911',  # copy_number_increase
    'part of contiguous gene duplication': 'SO_1000173'  # tandem_duplication
}

class G2P():

    def __init__(self):
        super(G2P, self).__init__()
        self.evidence_strings = list()
        self.unmapped_diseases = set()

        # Create OnToma object
        self.ontoma = ontoma.interface.OnToma()

    def map_disease_name_to_ontology(self, disease_name, omim_id):
        '''
        Searches the G2P disease name in EFO, ORDO, HP and MONDO

        OnToma is used to query the ontology OBO files, the manual mapping file and the Zooma and OLS APIs. MONDO is only searched if no exact matches are found.

        Args:
            disease_name (str): Gene2Phenotype disease name extracted from the "disease name" column
            omim_id (int): Gene2Phenotype disease OMIM id extracted from the "disease mim" column

        Returns:
            dict: Dictionary that contains id and name of mapped ontology term or `None` if not found
        '''

        logging.info(f"Mapping '{disease_name}'")

        # Search disease name using OnToma and accept perfect matches
        ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
        if ontoma_mapping:
            if ontoma_mapping['action'] is None:
                return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
            elif ontoma_mapping['quality'] == "match":
                # Match in HP or ORDO, check if there is a match in MONDO too. If so, give preference to MONDO hit
                logging.info(f"OnToma found match for '{disease_name}' in HP or ORDO, checking in MONDO")
                mondo_mapping = self.search_mondo(disease_name)
                if mondo_mapping:
                    if mondo_mapping['exact']:
                        logging.info(f"Using MONDO match")
                        return mondo_mapping
                    else:
                        logging.info(f"No exact matches in MONDO, using OnToma results")
                        return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
                else:
                    logging.info(f"No match in MONDO, using OnToma results")
                    return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
            else:
                # OnToma fuzzy match. First check if the mapping term has a xref to the OMIM id. If not, check in MONDO and if there is not match ignore evidence and report disease
                logging.info(f"Fuzzy match for '{disease_name}' found in OnToma, checking if it has a xref to OMIM:{omim_id}")
                if self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                    for efo_xref in self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                        # Extract EFO id from OnToma results
                        efo_id = ontoma_mapping['term'].split('/')[-1].replace('_', ':')
                        if efo_id == efo_xref['id']:
                            logging.info(
                                f"{ontoma_mapping['term']} has a xref to OMIM:{omim_id}, using this term as a match for {disease_name}")
                            return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
                # xref search didn't work, try MONDO as the last resort
                logging.info(f"No exact match for '{disease_name}' found in OnToma, checking in MONDO")
                mondo_mapping = self.search_mondo(disease_name)
                if mondo_mapping:
                    if mondo_mapping['exact']:
                        logging.info(f"Using MONDO match")
                        return mondo_mapping
                    else:
                        return
                else:
                    logging.info(f"Fuzzy match from OnToma ignored for '{disease_name}', recording the unmapped disease")
                    # Record the unmapped disease
                    self.unmapped_diseases.add((disease_name, True, ontoma_mapping['term'], ontoma_mapping['label']))
                    return
        else:
            # No match in EFO, HP or ORDO
            logging.info(f"No match found in EFO, HP or ORDO for '{disease_name}' with OnToma, checking in MONDO")
            mondo_mapping = self.search_mondo(disease_name)
            if mondo_mapping:
                if mondo_mapping['exact']:
                    logging.info(f"Using MONDO match")
                    return mondo_mapping
                else:
                    return
            else:
                logging.info(
                    f"'{disease_name}' could not be mapped to any EFO, HP, ORDO or MONDO. Skipping it, it should be checked with Gene2Phenotype and/or EFO")
                # Record the unmapped disease
                self.unmapped_diseases.add((disease_name, False, "", ""))
                return

    def search_mondo(self, disease_name):
        '''
        Searches the G2P disease name in MONDO, both using the OBO file and OLS

        Args:
            disease_name (str): Gene2Phenotype disease name extracted from the "disease name" column

        Returns:
            dict: Dictionary containing MONDO id, name and whether it's exact or not, or `None` if not found
        '''

        disease_name = disease_name.lower()

        # mondo_lookup works like a dictionary lookup so if disease is not in there it raises and error instead of returning `None`
        try:
            mondo_term = self.ontoma.mondo_lookup(disease_name)
            logging.info(f"Found {mondo_term} for '{disease_name}' from MONDO OBO file - Request EFO to import it from MONDO")
            return {'id': mondo_term, 'name': self.ontoma.get_mondo_label(mondo_term), 'exact': True}
        except KeyError as e:
            logging.info(f"No match found in MONDO OBO for '{disease_name}'")
            if self.ontoma._ols.besthit(disease_name, ontology=['mondo'], field_list=['iri', 'label'], exact=True):
                # Exact match
                exact_ols_mondo = self.ontoma._ols.besthit(disease_name, ontology=['mondo'],
                                                           field_list=['iri', 'label'],
                                                           exact=True)
                logging.info(
                    f"Found {exact_ols_mondo['iri']} as an exact match for '{disease_name}' from OLS API MONDO lookup - Request EFO to import it from MONDO")
                return {'id': exact_ols_mondo['iri'], 'name': exact_ols_mondo['label'], 'exact': True}
            elif self.ontoma._ols.besthit(disease_name, ontology=['mondo'], field_list=['iri', 'label'], bytype='class'):
                # Non-exact match, treat it as fuzzy matches above, i.e. ignore and report
                ols_mondo = self.ontoma._ols.besthit(disease_name,
                                                     ontology=['mondo'],
                                                     field_list=['iri', 'label'],
                                                     bytype='class')
                logging.info(f"Non-exact match from OnToma ignored for '{disease_name}', recording the unmapped disease")
                # Record the unmapped disease
                self.unmapped_diseases.add((disease_name, True, ols_mondo['iri'], ols_mondo['label']))
                return {'id': ols_mondo['iri'], 'name': ols_mondo['label'], 'exact': False}
            else:
                logging.info(f"No match for '{disease_name}' in MONDO")
                return

    def process_g2p(self, dd_file, eye_file, skin_file, cancer_file, evidence_file):

        # Parser DD file
        logging.info("Started parsing DD file")
        self.generate_evidence_strings(dd_file)
        # Parser eye file
        logging.info("Started parsing eye file")
        self.generate_evidence_strings(eye_file)
        # Parser skin file
        logging.info("Started parsing skin file")
        self.generate_evidence_strings(skin_file)
        # Parser cancer file
        logging.info("Started parsing cancer file")
        self.generate_evidence_strings(cancer_file)

        # Save results to file
        self.write_evidence_strings(evidence_file)

    def generate_evidence_strings(self, filename):

        total_efo = 0

        with gzip.open(filename, mode='rt') as zf:

            reader = csv.DictReader(zf, delimiter=',', quotechar='"')
            c = 0
            for row in reader:
                c += 1

                if c == 20:
                    break

                logging.info(f"Parsing row {c}")

                # Column names are:
                # "gene symbol","gene mim","disease name","disease mim","DDD category","allelic requirement",
                # "mutation consequence",phenotypes,"organ specificity list",pmids,panel,"prev symbols","hgnc id",
                # "gene disease pair entry date"
                gene_symbol = row["gene symbol"]
                disease_name = row["disease name"]
                disease_mim = row["disease mim"]
                allelic_requirement = row["allelic requirement"]
                mutation_consequence = row["mutation consequence"]
                confidence = row["DDD category"]
                pmids = row["pmids"]
                panel = row["panel"]

                # Map disease to ontology terms
                disease_mapping = self.map_disease_name_to_ontology(disease_name, disease_mim)
                if disease_mapping:
                    total_efo += 1

                    logging.info(f"{gene_symbol}, {gene_symbol}, '{disease_name}', {disease_mapping['id']}")

                    evidence = {
                        'datasourceId': 'gene2phenotype',
                        'datatypeId': 'genetic_literature',
                        'targetFromSourceId': gene_symbol.rstrip(),
                        'diseaseFromSource': disease_name,
                        'diseaseFromSourceId': disease_mim,
                        'diseaseFromSourceMappedId': ontoma.interface.make_uri(disease_mapping['id']).split("/")[-1],
                        'confidence': confidence,
                        'studyId': panel
                    }

                    # Add literature provenance if there are PMIDs
                    if len(pmids) != 0:
                        evidence['literature'] = pmids.split(";")

                    # Ignore allelic requirement if it's empty string
                    evidence['allelicRequirements'] = [allelic_requirement] if allelic_requirement else logging.warn('Empty allelic requirement, ignoring the value')

                    # Assign SO code based on mutation consequence field
                    if mutation_consequence in G2P_mutationCsq2functionalCsq:
                        if G2P_mutationCsq2functionalCsq[mutation_consequence]:
                            evidence['variantFunctionalConsequenceId'] = G2P_mutationCsq2functionalCsq[mutation_consequence]
                        else:
                            logging.error('Ignoring empty mutation consequence: {}'.format(mutation_consequence))
                    else:
                        logging.error('{} is not a recognised G2P mutation consequence, ignoring the value'.format(mutation_consequence))

                    self.evidence_strings.append(evidence)

            logging.info(f"Processed {c} diseases, mapped {total_efo}\n")

    def write_evidence_strings(self, filename):
        logging.info("Writing Gene2Phenotype evidence strings to %s", filename)
        with gzip.open(filename, 'w') as tp_file:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                json.dump(evidence_string, tp_file)
                tp_file.write("\n")
            logging.info(f"{n} evidence strings saved.\n")


def main(dd_file, eye_file, skin_file, cancer_file, outfile):

    g2p = G2P()
    g2p.process_g2p(dd_file, eye_file, skin_file, cancer_file, outfile)


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-d', '--dd_panel',
                        help='DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-e', '--eye_panel',
                        help='Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-k', '--skin_panel',
                        help='Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-c', '--cancer_panel',
                        help='Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-o', '--output_file', help='Name of gzipped evidence file', type=str)
    parser.add_argument('-l', '--log_file', help='Name of gzipped evidence file', type=str)

    args = parser.parse_args()

    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    outfile = args.output_file
    log_file = args.log_file

    # Configure logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Calling main:
    main(dd_file, eye_file, skin_file, cancer_file, outfile)
