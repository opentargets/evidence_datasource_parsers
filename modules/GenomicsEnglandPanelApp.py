from settings import Config
from common.HGNCParser import GeneParser
import python_jsonschema_objects as pjo
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from ontoma import OnToma
from datetime import datetime
import logging
import requests
import json
import argparse
import re


PanelApp_classification2score = {
    "green": 1,
    "amber": 0.5,
}

class PanelApp_evidence_generator():

    def __init__(self, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):

        # Get JSON Schema

        self.schema_version = schema_version
        schema_url = f"https://raw.githubusercontent.com/opentargets/json_schema/{self.schema_version}/draft4_schemas/opentargets.json"
        logging.info(f"Loading JSON Schema from {schema_url}")


        # Initialize JSON Schema builder

        try:
            r = requests.get(schema_url)
            r.raise_for_status()
            json_schema = r.json()
            self.builder = pjo.ObjectBuilder(json_schema)
            self.evidence_builder = self.builder.build_classes()
        except requests.exceptions.HTTPError as e:
            logging.error("Invalid JSON schema version")
            raise e

        # Create OnToma object

        self.otmap = OnToma()

        # Parse Gene Symbol to EnsEMBL ID

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

    @staticmethod
    def build_publications(dataframe):
        '''
        Builds a dataframe with the publications fetched from the PanelApp API and cleans them to match PubMed IDs.
        '''
        panel_lst = dataframe["Panel Id"].unique()
        lst = []

        for panel in panel_lst:
            tmp = dataframe[dataframe["Panel Id"] == panel]
            request = PanelApp_evidence_generator.publications_from_panel(panel)
            tmp["Publications"] = tmp.apply(lambda X: PanelApp_evidence_generator.publication_from_symbol(X.Symbol, request), axis=1)
            lst.append(tmp)
        
        dataframe = pd.concat(lst, ignore_index=True, sort=False)

        cleaned_publication = []
        for row in dataframe["Publications"].to_list():
            try:
                tmp = [re.match(r"(\d{8})", e)[0] for e in row]
                tmp = list(set(tmp)) # Removing duplicated publications
                cleaned_publication.append(tmp)
            except Exception as e:
                cleaned_publication.append([])
                continue

        dataframe["Publications"] = cleaned_publication

        return dataframe

    @staticmethod 
    def publications_from_panel(panel_id):
        '''
        Queries the PanelApp API to obtain a list of the publications for every gene within a panel_id
        '''
        try:
            url = f"http://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/"
            res = requests.get(url).json()

            return res["genes"]
        except:
            logging.error("Query of the PanelApp API failed.")
            return None
    
    @staticmethod
    def publication_from_symbol(symbol, response):
        '''
        Returns the list of publications for a given symbol in a PanelApp query response.
        '''
        for gene in response:
            if gene["gene_data"]["gene_symbol"] == symbol:
                publication = gene["publications"]
                return publication
            continue

    @staticmethod
    def clean_dataframe(dataframe):
        '''
        Cleaning of the initial .tsv or .csv file
        '''

        # NaNs and "No OMIM phenotype" in "Phenotypes" column --> Assignment of Pannel Name

        dataframe.loc[(dataframe["Phenotypes"] == "No OMIM phenotype") | (dataframe["Phenotypes"].isna()), "Phenotypes"] = dataframe["Panel Name"]

        # Handling multiple phenotypes column: Phenotypes exploded and phenotype_list with the original string

        dataframe["phenotype_list"] = dataframe["Phenotypes"]
        dataframe["Phenotypes"] = dataframe["Phenotypes"].apply(lambda X: X.split(";"))
        dataframe = dataframe.explode("Phenotypes")

        # Extracting and cleaning the OMIM codes: 
        #   removal of the OMIM codes in the Phenotypes column and the inclusion of these in a Codes column
        #   deleting special characters

        dataframe["Codes"] = dataframe["Phenotypes"].str.extract(r"(\d{6})")
        dataframe["Phenotypes"] = dataframe["Phenotypes"].replace("(\d{6})","",regex=True)

        dataframe["Phenotypes"] = dataframe["Phenotypes"].str.replace("-", " ", regex=False)
        dataframe["Phenotypes"] = dataframe["Phenotypes"].str.replace("[^0-9a-zA-Z *]", "", regex=True)
        dataframe["Phenotypes"] = dataframe["Phenotypes"].apply(lambda X: X.strip())

        return dataframe

    def ontoma_query(self, iterable, dict_name="ontoma_queries.json"):
        '''
        Queries the OnToma utility to map a phenotype to a disease.
        '''

        dct = dict()

        for e in iterable:
            try:
                tmp = self.otmap.find_term(e, verbose=True)
                if tmp != None:
                    dct[e] = tmp
                else:
                    dct[e] = {
                    'term': None,
                     'label': None,
                     'source': None,
                     'quality': None,
                     'action': None
                    }
            except Exception as error:
                logging.error(f"{e} mapping has failed.")
                continue
        
        with open(dict_name, "w") as outfile:  
            json.dump(dct, outfile)

        return dct

    @staticmethod
    def check_xref(phenotype, code, mappings_dict, codes_dict):
        '''
        Among the Fuzzy results of a phenotype query, it checks if the term and the respective code points to the same term
        '''

        try:
            if phenotype == "":
                return
            phenotype_result = mappings_dict[phenotype]
            if phenotype_result == None:
                return

            if phenotype_result["quality"] == "fuzzy":
                code_result = codes_dict[code]
                if code_result == None:
                    return

                if code_result["term"] == phenotype_result["term"]:
                    mappings_dict[phenotype]["quality"] = "match"
                    mappings_dict[phenotype]["action"] = "checked"
        except Exception as e:
            logging.error(f'No OMIM code for phenotype: {phenotype}')

    @staticmethod
    def build_mappings(mappings_dict, dataframe):
        '''
        Feeds the dataframe with the mappings results.
        '''
        # Splitting dictionaries between fuzzies and matches

        fuzzy = {}
        match = {}

        for i in mappings_dict.keys():
            if mappings_dict[i]["quality"] == "fuzzy":
                fuzzy[i] = mappings_dict[i]
            elif mappings_dict[i]["quality"] == "match":
                match[i] = mappings_dict[i]

        # Creating the corresponding OnToma result for each phenotype
        # New columns: OnToma Result, OnToma Source, OnToma Term, OnToma Label

        for phenotype in match.keys():
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Result"] = "match"
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Action"] = match[phenotype]["action"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Term"] = match[phenotype]["term"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Label"] = match[phenotype]["label"]
            
        for phenotype in fuzzy.keys():
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Result"] = "fuzzy"
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Action"] = fuzzy[phenotype]["action"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Term"] = fuzzy[phenotype]["term"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Label"] = fuzzy[phenotype]["label"]

        return dataframe

    def get_pub_array(self):
        '''
        Returns an array where all the publications are listed as dictionaries
        '''
        pub_array = []

        if len(self.publication) != 0:
            for pub in self.publication:
                pub_dict = {
                    "lit_id": f"http://europepmc.org/abstract/MED/{pub}"
                }
                pub_array.append(pub_dict)
        else:
            pub_dict = {
                "lit_id": None
            }
            pub_array.append(pub_dict)
        return pub_array
    

    def get_evidence_string(self, row):


        # Association level information:
        self.gene_symbol = row["Symbol"]
        self.mapped_disease = row["OnToma Label"]
        self.mapped_id = row["OnToma Term"]
        self.source_array = row["phenotype_list"].split(";")
        self.source_name = row["Phenotypes"]
        self.mode_of_inheritance = row["Mode of inheritance"]
        self.evidence_classification = row["List"]
        self.sources = row["Sources"]
        self.publication = row["Publications"]
        
        try:
            self.ensembl_iri = "http://identifiers.org/ensembl/" + self.genes[self.gene_symbol]
        except:
            self.ensembl_iri = "http://identifiers.org/ensembl/" + row["EnsemblId(GRch37)"]

        # Panel level information:
        self.panel_id = row["Panel Id"]
        self.panel_name = row["Panel Name"]
        self.panel_version = row["Panel Version"]
        self.panel_type = row["Panel Types"]
        self.super_panel_id = row["Super Panel Id"]
        self.super_panel_name = row["Super Panel Name"]
        self.super_panel_version = row["Super Panel Version"]

        # If the association has no EFO attached:
        if self.mapped_id is None or self.mapped_id == '':
            logging.warning(f'No EFO id for association row: {row.name}')
            return None

        target_field = {
                'id' : self.ensembl_iri,
                'activity' : "http://identifiers.org/cttv.activity/unknown",
                'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                'target_name' : self.gene_symbol
        }

        provenance_type = {
                        'database' : {
                            'id' : "Genomics England PanelApp",
                            'version' : datetime.now().isoformat(),
                            'dbxref' : {
                                'id' : "Genomics England PanelApp",
                                'url': "https://panelapp.genomicsengland.co.uk/",
                                'version' : datetime.now().isoformat() # Config.GE_PANEL_VERSION
                                        }
                                    }
                        }
        # literature field only shown when the publication array is not empty
        if len(self.publication) != 0:
            provenance_type["literature"] = {
                            'references': self.get_pub_array()
                        }

        resource_score = {
                        'method': {
                            'description': "Further details in the Genomics England PanelApp.",
                            'url': 'https://panelapp.genomicsengland.co.uk'
                            },
                        'type': "probability",
                        'value': PanelApp_classification2score[self.evidence_classification]
                    }
        
        urls = [
                {
                'nice_name' : f'Further details in the Genomics England PanelApp for panel {self.panel_id} and gene {self.gene_symbol}',
                'url' : f"https://panelapp.genomicsengland.co.uk/panels/{self.panel_id}/{self.gene_symbol}"
                }
            ]

        evidence_field = {
                    'date_asserted' : datetime.now().isoformat(),
                    'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000205"],
                    'is_associated' : True,
                    'provenance_type' : provenance_type,
                    'resource_score' : resource_score,
                    'urls' : urls,
                    'confidence' : self.evidence_classification,
                    'panel_id': self.panel_id,
                    'panel_version': self.panel_version,
                    'allelic_requirement' : self.mode_of_inheritance,
                    'phenotypes_array' : self.phenotypes_array
                    }

        target_field = {
                        'id' : self.ensembl_iri,
                        'activity' : "http://identifiers.org/cttv.activity/predicted_damaging",
                        'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                        'target_name' : self.gene_symbol
                    }

        disease_field = {
                        'id': self.mapped_id,
                        'name' : self.mapped_disease,
                        'source_name': self.source_name
                    }
        
        unique_association_field = {
                        'disease_iri': self.mapped_id,
                        'target_id': self.ensembl_iri,
                        'panel_id': self.panel_id,
                        'panel_name': self.panel_name,
                        'panel_version': self.panel_version,
                        'original_disease_name': self.source_disease,
                    }

        try:
            evidence = self.evidence_builder.Opentargets5(
                type = "genetic_literature",
                access_level = "public",
                sourceID = "genomics_england",
                evidence = evidence_field,
                target = target_field,
                disease = disease_field,
                unique_association_fields = unique_association_field,
                validated_against_schema_version = self.schema_version
            )
            return evidence.serialize()
        except Exception as e:
            print(e)
            logging.error(f'Evidence generation failed for row: {row.name}')
            raise

    def write_evidence_strings(self, dataframe, mappings_dict):

        # Read input file
        dataframe = pd.read_csv(dataframe, sep='\t')

        # Filtering with green and ambar evidences

        dataframe = dataframe[((dataframe["List"] == "green") | (dataframe["List"] == "amber")) & (dataframe["Panel Version"] > 1) & (dataframe["Panel Status"] == "PUBLIC")].reset_index(drop=True)
        dataframe.dropna(axis=1, how="all", inplace=True)

        # Feed the dataframe with publications

        logging.info("Fetching publications from the API...")
        dataframe = self.build_publications(dataframe)
        logging.info("Publications loaded.")

        # Splitting and cleaning the dataframe according to the phenotype string
        
        dataframe = PanelApp_evidence_generator.clean_dataframe(dataframe)

        # Mapping the phenotypes and OMIM codes to an EFO code

        codes = dataframe["Codes"].unique()
        phenotypes = dataframe["Phenotypes"].unique()

        if len(mappings_dict) == 0:
            mappings_dict = self.ontoma_query(phenotypes)
            logging.info("Disease mappings completed.")
        else:
            logging.info("Disease mappings imported.")
        
        codes_dict = self.ontoma_query(codes)

        # Xref of the fuzzy matches with the corresponding code

        phenotypes_list = dataframe["Phenotypes"].to_list()
        codes_list = dataframe["Codes"].to_list()
        merged = list(zip(phenotypes_list, codes_list))

        for pheno, code in merged:
            PanelApp_evidence_generator.check_xref(pheno, code, mappings_dict, codes_dict)

        # New columns: OnToma Result, OnToma Source, OnToma Term, OnToma Label
        dataframe = PanelApp_evidence_generator.build_mappings(mappings_dict, dataframe)
        dataframe = dataframe[dataframe["OnToma Result"] == "match"]

        # Build evidence strings per row
        evidences = dataframe.apply(self.get_evidence_string, axis=1)
        logging.info(f"{len(evidences)} evidence strings have been generated.")

        return evidences

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates Genomics England PanelApp sourced evidences.")

    parser.add_argument("-i", "--input_file", required=True, type=str, help="Input .tsv file with the table containing association details.")
    parser.add_argument("-o", "--output_file", required=True, type=str, help="Name of the json output file containing the evidence strings.")
    parser.add_argument("-s", "--schema_version", required=True, type=str, help="JSON schema version to use, e.g. 1.6.9. It must be branch or a tag available in https://github.com/opentargets/json_schema.")
    parser.add_argument("-d", "--dictionary", required=False, type=str, help="Path of the dictionary containing the mapped phenotypes.")

    # Parsing parameters
    args = parser.parse_args()

    dataframe = args.input_file
    output_file = args.output_file
    schema_version = args.schema_version

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidence_builder = PanelApp_evidence_generator(schema_version)

    # Import dictionary if present
    if args.dictionary:
        with open(args.dictionary, "r") as f:
                mappings_dict = json.load(f)
    else:
        mappings_dict = {}
    
    # Writing evidence strings into a json file
    evidences = evidence_builder.write_evidence_strings(dataframe, mappings_dict)

    with open(output_file, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n'))
    
    logging.info(f"Evidence strings saved into {output_file}. Exiting.")

if __name__ == '__main__':
    main()