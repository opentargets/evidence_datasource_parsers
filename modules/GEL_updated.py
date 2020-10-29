import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from ontoma import OnToma

from settings import Config
from common.HGNCParser import GeneParser
from datetime import datetime
import logging
import requests
import json
import python_jsonschema_objects as pjo
import argparse


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
            logging.error('Invalid JSON schema version')
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
        Builds a dataframe with the publications fetched from the PanelApp API
        '''
        panel_lst = dataframe["Panel Id"].unique()
        lst = []

        for panel in panel_lst:
            tmp = dataframe[dataframe["Panel Id"] == panel]
            request = publications_from_panel(panel)
            tmp["Publications"] = tmp.apply(lambda X: publication_from_symbol(X.Symbol, request), axis=1)
            lst.append(tmp)
        
        dataframe = pd.concat(lst, ignore_index=True, sort=False)
        return dataframe

    @staticmethod 
    def publications_from_panel(panel_id):
        '''
        queries the PanelApp API to obtain a list of the publications for every gene within a panel_id
        '''
        try:
            url = f"http://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/"
            res = requests.get(url).json()

            return res["genes"]
        except:
            return None
    
    @staticmethod
    def publication_from_symbol(symbol, response):
        '''

        '''
        for gene in response:
            if gene["gene_data"]["gene_symbol"] == symbol:
                publication = gene["publications"]
                return publication
            continue

    @staticmethod
    def split_dataframes(dataframe):
        '''
        Cleaning of the initial .tsv or .csv file
        '''

        # NaNs in "Phenotypes" column --> Assignment of Pannel Name

        dataframe.loc[dataframe["Phenotypes"].isna(), "Phenotypes"] = dataframe["Panel Name"]

        # Splitting between one and multiple associations 
        # "Phenotypes" column contains the separator ";"

        one_phenotype = dataframe[~dataframe["Phenotypes"].str.contains(";", na=False)].reset_index(drop=True)
        multiple_phenotype = dataframe[dataframe["Phenotypes"].str.contains(";", na=False)].reset_index(drop=True)

        # Among one_phenotype, splitting between those that include an OMIM code and those who don't.

        OMIM_codes = one_phenotype[one_phenotype["Phenotypes"].str.contains('[0-9]{5}', na=False, regex=True)].reset_index(drop=True)
        not_OMIM = one_phenotype[~one_phenotype["Phenotypes"].str.contains('[0-9]{5}', na=False, regex=True)].reset_index(drop=True)

        # Cleaning the OMIM_codes df: 
        #   - cleaned_OMIM: df after the removal of the OMIM codes in the Phenotypes column and the inclusion of these in a Codes column
        #   - Deleting special characters (Code is crashing, solution is a shame - TO FIX)

        cleaned_OMIM = OMIM_codes.copy()
        cleaned_OMIM["Codes"] = cleaned_OMIM["Phenotypes"].str.extract(r"(\d{6})")
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].replace("(\d{6})","",regex=True)

        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("{", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("}", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(",", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("(", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(")", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("[", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("]", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("'", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace('"', ''))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(".", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("-", " "))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("?", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("#", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.strip())

        # Exporting dfs to .csv

        #cleaned_OMIM.to_csv("cleaned_OMIM.csv", index=False)
        #not_OMIM.to_csv("not_OMIM.csv", index=False)
        #multiple_phenotype.to_csv("multiple_phenotype.csv", index=False)

        return cleaned_OMIM, not_OMIM, multiple_phenotype

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
                print(error)
                continue
        
        with open(dict_name, "w") as outfile:  
            json.dump(dct, outfile)

        return dct
    
    @staticmethod
    def build_mappings(mappings_dict, dataframe):

        '''valid_mappings = {}
        for i in mappings_dict.keys():
            if mappings_dict[i]["quality"] != None:
                valid_mappings[i] = mappings_dict[i]

        for phenotype in valid_mappings.keys():
            try:
                data.loc[data["Phenotypes"] == phenotype, "OnToma Result"] = valid_mappings[phenotype]["quality"]
                data.loc[data["Phenotypes"] == phenotype, "OnToma Source"] = valid_mappings[phenotype]["source"]
                data.loc[data["Phenotypes"] == phenotype, "OnToma Term"] = valid_mappings[phenotype]["term"]
                data.loc[data["Phenotypes"] == phenotype, "OnToma Label"] = valid_mappings[phenotype]["label"]
            except:
                continue'''
        
        # Splitting dictionaries between fuzzies and matches

        fuzzy = {}
        match = {}

        for i in mappings_dict.keys():
            try:
                if mappings_dict[i]["quality"] == "fuzzy":
                    fuzzy[i] = mappings_dict[i]
                elif mappings_dict[i]["quality"] == "match":
                    match[i] = mappings_dict[i]
            except:
                continue

        # Creating the corresponding OnToma result for each phenotype
        # New columns: OnToma Result, OnToma Source, OnToma Term, OnToma Label
        # TO-DO: Refactor this

        for phenotype in match.keys():
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Result"] = "match"
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Source"] = match[phenotype]["source"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Term"] = match[phenotype]["term"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Label"] = match[phenotype]["label"]
            
        for phenotype in fuzzy.keys():
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Result"] = "fuzzy"
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Source"] = fuzzy[phenotype]["source"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Term"] = fuzzy[phenotype]["term"]
            dataframe.loc[dataframe["Phenotypes"] == phenotype, "OnToma Label"] = fuzzy[phenotype]["label"]

        return dataframe
    

    def get_evidence_string(self, row):


        # Association level information:
        self.gene_symbol = row["Symbol"]
        self.ensembl_iri = "http://identifiers.org/ensembl/" + self.genes[self.gene_symbol]
        self.mapped_disease = row["OnToma Label"]
        self.mapped_id = row["OnToma Term"]
        self.source_disease = row["Phenotypes"]
        self.mode_of_inheritance = row["Mode of inheritance"]
        self.evidence_classification = row["List"]
        self.sources = row["Sources"]
        #self.publication = row["Publications"]

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
                    'urls' : urls
                    #'confidence' : PanelApp_classification2score[self.evidence_classification]
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
                        'source_name': self.source_disease
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
                type = "genetic_literature", # done
                access_level = "public", # done
                sourceID = "genomics_england", # done
                evidence = evidence_field, # done
                target = target_field, # done
                disease = disease_field, # done
                unique_association_fields = unique_association_field, # done
                validated_against_schema_version = self.schema_version
            )
            return evidence.serialize()
        except:
            logging.warning(f'Evidence generation failed for row: {row.name}')
            raise

def write_evidence_strings(dataframe, schema_version):

    # Initialize logging:
    logging.basicConfig(
        filename='evidence_builder.log',
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidence_builder = PanelApp_evidence_generator(schema_version)

    # Read input file
    dataframe = pd.read_csv(dataframe, sep='\t')

    # Filtering with green and ambar evidences

    dataframe = dataframe[((dataframe["List"] == "green") | (dataframe["List"] == "amber")) & (dataframe["Panel Version"] > 1) & (dataframe["Panel Status"] == "PUBLIC")].reset_index(drop=True)
    dataframe.dropna(axis=1, how="all", inplace=True)

    # Feed the dataframe with publications

    # dataframe = evidence_builder.build_publications(dataframe)

    # Splitting and cleaning the dataframe according to the phenotype string

    OMIM_codes, not_OMIM, multiple_phenotype = evidence_builder.split_dataframes(dataframe)

    # Mapping the phenotypes to an EFO code (Excluding multiple phenotypes ATM)

    dataframe = pd.concat([OMIM_codes, not_OMIM], ignore_index=True)

    phenotypes = dataframe["Phenotypes"].unique()

    mappings_dict = evidence_builder.ontoma_query(phenotypes)
    assert mappings_dict, "Disease mappings completed."

    # New columns: OnToma Result, OnToma Source, OnToma Term, OnToma Label
    dataframe = evidence_builder.build_mappings(mappings_dict,dataframe)

    # Build evidence strings per row
    evidences = dataframe.apply(evidence_builder.get_evidence_string, axis=1)

    return evidences

def main():
    '''

    '''
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates Genomics England PanelApp sourced evidences.")

    parser.add_argument("--input_file", required=True, type=str, help="Input .tsv file with the table containing association details.")
    parser.add_argument("--output_file", required=True, type=str, help="Name of the json output file containing the evidence strings.")
    parser.add_argument("--schema_version", required=True, type=str, help="JSON schema version to use, e.g. 1.6.9. It must be branch or a tag available in https://github.com/opentargets/json_schema.")

    # Parsing parameters
    args = parser.parse_args()

    dataframe = args.input_file
    output_file = args.output_file
    schema_version = args.schema_version

    # Writing evidence strings into a json file
    evidences = write_evidence_strings(dataframe, schema_version)

    with open(output_file, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n'))

if __name__ == '__main__':
    
    main()