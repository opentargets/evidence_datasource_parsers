from settings import Config
import ontoma

import python_jsonschema_objects as pjo

import logging
import csv
import gzip
import requests
import argparse
import json

G2P_mutationCsq2functionalCsq = {
    'loss of function' : 'SO_0002054', #loss_of_function_variant
    'all missense/in frame' : 'SO_0001650', #inframe_variant
    'uncertain' : 'SO_0002220', #function_uncertain_variant
    'activating' : 'SO_0002053', #gain_of_function_variant
    'dominant negative' : 'SO_0002052', #dominant_negative_variant
    '' : None,
    'gain of function' : 'SO_0002053', #gain_of_function_variant
    'cis-regulatory or promotor mutation' : 'SO_0001566', #regulatory_region_variant
    '5_prime or 3_prime UTR mutation' : 'SO_0001622', #UTR_variant
    'increased gene dosage' : 'SO_0001911', #copy_number_increase
    'part of contiguous gene duplication' : 'SO_1000173' #tandem_duplication
}



class G2P():
    def __init__(self, g2p_version):
        super(G2P, self).__init__()
        self.evidence_strings = list()
        self.unmapped_diseases = set()
        self.g2p_version = g2p_version

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

        self._logger.info(f"Mapping '{disease_name}'")

        # Search disease name using OnToma and accept perfect matches
        ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
        if ontoma_mapping:
            if ontoma_mapping['action'] is None:
                return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
            elif ontoma_mapping['quality'] == "match":
                # Match in HP or ORDO, check if there is a match in MONDO too. If so, give preference to MONDO hit
                self._logger.info(f"OnToma found match for '{disease_name}' in HP or ORDO, checking in MONDO")
                mondo_mapping = self.search_mondo(disease_name)
                if mondo_mapping:
                    if mondo_mapping['exact']:
                        self._logger.info(f"Using MONDO match")
                        return mondo_mapping
                    else:
                        self._logger.info(f"No exact matches in MONDO, using OnToma results")
                        return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
                else:
                    self._logger.info(f"No match in MONDO, using OnToma results")
                    return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
            else:
                # OnToma fuzzy match. First check if the mapping term has a xref to the OMIM id. If not, check in MONDO and if there is not match ignore evidence and report disease
                self._logger.info(f"Fuzzy match for '{disease_name}' found in OnToma, checking if it has a xref to OMIM:{omim_id}")
                if self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                    for efo_xref in self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                        # Extract EFO id from OnToma results
                        efo_id = ontoma_mapping['term'].split('/')[-1].replace('_', ':')
                        if efo_id == efo_xref['id']:
                            self._logger.info(
                                f"{ontoma_mapping['term']} has a xref to OMIM:{omim_id}, using this term as a match for {disease_name}")
                            return {'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}
                # xref search didn't work, try MONDO as the last resort
                self._logger.info(f"No exact match for '{disease_name}' found in OnToma, checking in MONDO")
                mondo_mapping = self.search_mondo(disease_name)
                if mondo_mapping:
                    if mondo_mapping['exact']:
                        self._logger.info(f"Using MONDO match")
                        return mondo_mapping
                    else:
                        return
                else:
                    self._logger.info(f"Fuzzy match from OnToma ignored for '{disease_name}', recording the unmapped disease")
                    # Record the unmapped disease
                    self.unmapped_diseases.add((disease_name, True, ontoma_mapping['term'], ontoma_mapping['label']))
                    return
        else:
            # No match in EFO, HP or ORDO
            self._logger.info(f"No match found in EFO, HP or ORDO for '{disease_name}' with OnToma, checking in MONDO")
            mondo_mapping = self.search_mondo(disease_name)
            if mondo_mapping:
                if mondo_mapping['exact']:
                    self._logger.info(f"Using MONDO match")
                    return mondo_mapping
                else:
                    return
            else:
                self._logger.info(
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
            self._logger.info(f"Found {mondo_term} for '{disease_name}' from MONDO OBO file - Request EFO to import it from MONDO")
            return {'id': mondo_term, 'name': self.ontoma.get_mondo_label(mondo_term), 'exact': True}
        except KeyError as e:
            self._logger.info(f"No match found in MONDO OBO for '{disease_name}'")
            if self.ontoma._ols.besthit(disease_name, ontology=['mondo'], field_list=['iri', 'label'], exact=True):
                # Exact match
                exact_ols_mondo = self.ontoma._ols.besthit(disease_name, ontology=['mondo'],
                                                           field_list=['iri', 'label'],
                                                           exact=True)
                self._logger.info(
                f"Found {exact_ols_mondo['iri']} as an exact match for '{disease_name}' from OLS API MONDO lookup - Request EFO to import it from MONDO")
                return {'id': exact_ols_mondo['iri'], 'name': exact_ols_mondo['label'], 'exact':True}
            elif self.ontoma._ols.besthit(disease_name, ontology=['mondo'], field_list=['iri', 'label'], bytype='class'):
                # Non-exact match, treat it as fuzzy matches above, i.e. ignore and report
                ols_mondo = self.ontoma._ols.besthit(disease_name,
                                                     ontology=['mondo'],
                                                     field_list=['iri', 'label'],
                                                     bytype='class')
                self._logger.info(f"Non-exact match from OnToma ignored for '{disease_name}', recording the unmapped disease")
                # Record the unmapped disease
                self.unmapped_diseases.add((disease_name, True, ols_mondo['iri'], ols_mondo['label']))
                return {'id': ols_mondo['iri'], 'name': ols_mondo['label'], 'exact': False}
            else:
                self._logger.info(f"No match for '{disease_name}' in MONDO")
                return

    def process_g2p(self, dd_file, eye_file, skin_file, cancer_file, evidence_file, unmapped_diseases_filename):

        # Parser DD file
        self._logger.info("Started parsing DD file")
        self.generate_evidence_strings(dd_file)
        # Parser eye file
        self._logger.info("Started parsing eye file")
        self.generate_evidence_strings(eye_file)
        # Parser skin file
        self._logger.info("Started parsing skin file")
        self.generate_evidence_strings(skin_file)
        # Parser cancer file
        self._logger.info("Started parsing cancer file")
        self.generate_evidence_strings(cancer_file)

        # Save results to file
        self.write_evidence_strings(evidence_file)

        # If selected, write unmapped diseases to file
        if unmapped_diseases_filename:
            self.write_unmapped_diseases(unmapped_diseases_filename)

    def generate_evidence_strings(self, filename):

        total_efo = 0

        with gzip.open(filename, mode='rt') as zf:

            reader = csv.DictReader(zf, delimiter=',', quotechar='"')
            c = 0
            for row in reader:
                c += 1
                self._logger.info(f"Parsing row {c}")

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
                    total_efo +=1

                    self._logger.info(f"{gene_symbol}, {target}, '{disease_name}', {disease_mapping['id']}")

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
                    evidence['allelicRequirements'] = [allelic_requirement] if allelic_requirement else self._logger.warn('Empty allelic requirement, ignoring the value')

                    # Assign SO code based on mutation consequence field
                    if mutation_consequence in G2P_mutationCsq2functionalCsq:
                        evidence['variantFunctionalConsequenceId'] = G2P_mutationCsq2functionalCsq[mutation_consequence]
                    else:
                        self._logger.error( '{} is not a recognised G2P mutation consequence, ignoring the value'.format(mutation_consequence))

                    self.evidence_strings.append(evidence)

            self._logger.info(f"Processed {c} diseases, mapped {total_efo}\n")

    def write_evidence_strings(self, filename):
        self._logger.info("Writing Gene2Phenotype evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                json.dump(evidence_string, tp_file)
                tp_file.write("\n")
            self._logger.info(f"{n} evidence strings saved.\n")

    def write_unmapped_diseases(self, filename):
        self._logger.info(f"Writing Gene2Phenotype diseases not mapped to EFO to {filename}")
        with open(filename, 'w') as unmapped_diseases_file:
            unmapped_diseases_file.write("disease_name\tontoma_mapping_bool\tontoma_disease_id\tontoma_disease_name\n")
            for unmapped_disease in self.unmapped_diseases:
                unmapped_diseases_file.write(unmapped_disease[0] + "\t" + str(unmapped_disease[1]) + "\t" + unmapped_disease[2] + "\t" + unmapped_disease[3] + "\n")


def main():

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-d', '--dd_panel',
                        help='DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str, default=Config.G2P_DD_FILENAME)
    parser.add_argument('-e', '--eye_panel',
                        help='Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str, default=Config.G2P_eye_FILENAME)
    parser.add_argument('-k', '--skin_panel',
                        help='Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str, default=Config.G2P_skin_FILENAME)
    parser.add_argument('-c', '--cancer_panel',
                        help='Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str, default=Config.G2P_cancer_FILENAME)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file. It uses the value of G2P_EVIDENCE_FILENAME in setting.py if not specified',
                        type=str, default=Config.G2P_EVIDENCE_FILENAME)
    parser.add_argument('-u', '--unmapped_diseases_file',
                        help='If specified, the diseases not mapped to EFO will be stored in this file',
                        type=str, default=False)
    parser.add_argument('-v', '--g2p_version',
                        help='Version of the Gene2Phenotype data used. If not available please use the date in which the data was downloaded in YYYY-MM-DD format',
                        type=str, required=True)

    args = parser.parse_args()
    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    outfile = args.output_file
    unmapped_diseases_file = args.unmapped_diseases_file
    g2p_version = args.g2p_version

    g2p = G2P(g2p_version=g2p_version)
    g2p.process_g2p(dd_file, eye_file, skin_file, cancer_file, outfile, unmapped_diseases_file)

if __name__ == "__main__":
    main()