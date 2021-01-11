from common.HGNCParser import GeneParser
from ontoma import OnToma
import logging
import requests
import argparse
import re
from multiprocessing import Pool
import functools as fn
import numpy as np
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

class PanelAppEvidenceGenerator():
    def __init__(self, dataframe, mappingStep, mappingsDict):
        # Initialize mapping variables
        self.mappingStep = mappingStep
        self.phenotypesMappings = mappingsDict
        self.codesMappings = None

        # Create OnToma object
        self.otmap = OnToma()

        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('evidence_builder') \
                .getOrCreate()

        self.dataframe = dataframe

    def writeEvidenceFromSource(self):
        '''
        Processing of the dataframe to build all the evidences from its data

        Args:
            dataframe (pyspark.DataFrame): Initial .tsv file
        Returns:
            evidences (array): Object with all the generated evidences strings
        '''

        # Read input file
        self.dataframe = self.spark.read.csv(self.dataframe, sep=r'\t', header=True)

        # Filtering with green and ambar evidences
        self.dataframe = self.dataframe.filter(
                                ((col("List") == "green" ) | (col("List") == "amber"))
                                                            &
                                (col("Panel Version") > 1)
                                                            &
                                (col("Panel Status") == "PUBLIC")                
        )

        # TODO: Feed the dataframe with publications
        logging.info("Fetching publications from the API...")
        #dataframe = self.buildPublications(dataframe)
        logging.info("Publications loaded.")

        # Splitting and cleaning the dataframe according to the phenotype string
        self.dataframe = PanelAppEvidenceGenerator.cleanDataframe(self.dataframe)

        # Map the diseases to an EFO term if necessary
        if self.mappingStep:
            self.dataframe = self.runMappingStep()

        # Build evidence strings per row
        evidences = dataframe.rdd.map(
            PanelAppEvidenceGenerator.parseEvidenceString)
            .collect() # list of dictionaries
        logging.info(f"{len(evidences)} evidence strings have been generated.")

        # WARNING! Given the explosion of phenotypes, it is necessary to remove redundant evidences
        #evidences = PanelAppEvidenceGenerator.removeRedundancies(evidences)

        return evidences

    def buildPublications(self):
        '''
        Populates a dataframe with the publications fetched from the PanelApp API and cleans them to match PubMed IDs.

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pyspark.DataFrame): DataFrame with an additional column: Publications
        '''
        pass

    @staticmethod 
    def publicationsFromPanel(panelId):
        '''
        Queries the PanelApp API to obtain a list of the publications for every gene within a panelId
        
        Args:
            panelId (str): Panel ID extracted from the "Panel Id" column
        Returns:
            response (dict): Response of the API containing all genes related to a panel and their publications
        '''
        try:
            url = f"http://panelapp.genomicsengland.co.uk/api/v1/panels/{panelId}/"
            res = requests.get(url).json()
            return res["genes"]
        except:
            logging.error("Query of the PanelApp API failed.")
            return None
    
    @staticmethod
    def publicationFromSymbol(symbol, response):
        '''
        Returns the list of publications for a given symbol in a PanelApp query response.
        Args:
            symbol (str): Gene symbol extracted from the "Symbol" column
            response (dict): Response of the API containing all genes related to a panel and their publications
        Returns:
            publication (list): Array with all publications for a particular gene in the corresponding Panel ID 
        '''
        for gene in response:
            if gene["gene_data"]["gene_symbol"] == symbol:
                publication = gene["publications"]
                return publication


    @staticmethod
    def cleanDataframe(dataframe):
        '''
        Args:
            dataframe (pyspark.DataFrame): Initial .tsv data
        Returns:
            dataframe (pyspark.DataFrame): Transformed initial dataframe
        '''
        # TODO: Wrap different steps into one function to iterate less times
        
        # NaNs and "No OMIM phenotype" in "Phenotypes" column --> Assignment of Pannel Name
        dataframe = dataframe.withColumn("cohortPhenotypes", when(col("Phenotypes") == "No OMIM phenotype", col("Panel Name"))
                                 .when(col("Phenotypes").isNull(), col("Panel Name"))
                                 .otherwise(col("Phenotypes")))
        
        # Handling multiple phenotypes column: 
        #   cohortPhenotypes --> original string 
        #   phenotype --> explosion of cohortPhenotypes
        splitLambda = udf(
            lambda X: X.split(";"),
            ArrayType(StringType()))
        dataframe = dataframe \
                            .withColumn("phenotype", splitLambda(col("cohortPhenotypes"))) \
                            .withColumn('phenotype', explode('phenotype'))

        # Extracting and cleaning the OMIM codes: 
        #   removal of the OMIM codes in the Phenotypes column and the inclusion in omimCodes
        #   deleting special characters
        stripLambda = udf(
            lambda X: X.strip(),
            StringType())
        dataframe = dataframe \
                        .withColumn("omimCode", regexp_extract(col("phenotype"), "(\d{6})", 1)) \
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "(\d{6})", "")) \
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "[^0-9a-zA-Z *]", "")) \
                        .withColumn("phenotype", stripLambda(col("phenotype")))

        return dataframe
    
    def runMappingStep(self):
        '''
        Builds the mapping logic into the dataframe

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pyspark.DataFrame): DataFrame with the mapping results filtered by only matches
        '''
        # Mapping the phenotypes and OMIM codes to an EFO code - 2 steps:
        #   Querying OnToma with all distinch codes and phenotypes
        #   Xref between the results of every phenotype/OMIM code pair for a better coverage
        omimCodesDistinct = list(self.dataframe.select("omimCode").distinct().toPandas()["omimCode"])
        phenotypesDistinct = list(self.dataframe.select("phenotype").distinct().toPandas()["omimCode"])
        omimCodes = list(self.dataframe.toPandas()["omimCode"])
        phenotypes = list(self.dataframe.toPandas()["phenotype"])
        phenotypeCodePairs = list(zip(phenotypes, omimCodes))

        if len(self.phenotypesMappings) == 0:
            # Checks whether the dictionary is not provided as a parameter
            try:
                pool = Pool(mp.cpu_count())
                self.phenotypesMappings = pool.map(self.diseaseToEfo, phenotypesDistinct) # list of dictionaries
                self.phenotypesMappings = {k:v for dct in self.phenotypesMappings for (k,v) in self.phenotypesMappings.items()}
                pool.close()
            except:
                self.phenotypesMappings = self.diseaseToEfo(*phenotypesDistinct)
            logging.info("Disease mappings completed.")
        else:
            logging.info("Disease mappings imported.")
        self.codesMappings = self.diseaseToEfo(*omimCodesDistinct) # TODO: add posibility to provide dict

        for pheno, code in phenotypeCodePairs:
            self.phenotypeCodePairCheck(pheno, code)

        # TODO: Add new columns: OnToma Result, OnToma Term, OnToma Label
        #self.dataframe = self.buildMappings()

        return self.dataframe.filter(col("ontomaResult") == "match")

    def diseaseToEfo(self, *iterable, dictExport="diseaseToEfo_results.json"):
        '''
        Queries the OnToma utility to map a phenotype to a disease.
        OnToma is used to query the ontology OBO files, the manual mapping file and the Zooma and OLS APIs.

        Args:
            iterable (array): Array or column of a dataframe containing the strings to query
            dictExport (str): Name of the output file where the OnToma queries will be saved
        Returns:
            mappings (dict): Output file. Keys: queried term (phenotype or OMIM code), Values: OnToma output
        '''
        mappings = dict()
        for e in iterable:
            try:
                ontomaResult = self.otmap.find_term(e, verbose=True)
                if ontomaResult != None:
                    mappings[e] = ontomaResult
                else:
                    mappings[e] = {
                    'term': None,
                     'label': None,
                     'source': None,
                     'quality': None,
                     'action': None
                    }
            except Exception as error:
                logging.error(f"{e} mapping has failed.")
                continue
        
        with open(dictExport, "w") as outfile:
            json.dump(mappings, outfile)
        
        return mappings

    def phenotypeCodePairCheck(self, phenotype, omimCode):
        '''
        Among the Fuzzy results of a phenotype query, it checks if the phenotype and the respective code points to the same EFO term

        Args:
            phenotype (str): Phenotype of each phenotype-OMIM code pair
            omimCode (str): Code of each phenotype-OMIM code pair
        '''
        try:
            if phenotype == "":
                return
            phenotypeResult = self.phenotypesMappings[phenotype]
            if phenotypeResult == None:
                return

            if phenotypeResult["quality"] == "fuzzy":
                codeResult = self.codesMappings[code]
                if codeResult == None:
                    return

                if codeResult["term"] == phenotypeResult["term"]:
                    # If both EFO terms coincide, the phenotype in mappingsDict becomes a match
                    self.phenotypesMappings[phenotype]["quality"] = "match"
                    self.phenotypesMappings[phenotype]["action"] = "checked"
        except Exception as e:
            logging.error(f'No OMIM code for phenotype: {phenotype}')

    def buildMappings(dataframe):
        '''
        Populates the dataframe with the mappings resulted from OnToma.

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
        Return:
            dataframe (pyspark.DataFrame): DataFrame with new columns corresponding to the OnToma result
        '''

        pass
    
    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "genomics_england",
                "datatypeId" : "genetic_literature",
                "confidence" : row["List"],
                "diseaseFromSource" : row["phenotype"],
                "diseaseFromSourceId" : row["omimCode"],
                "diseaseFromSourceMappedId" : row["ontomaId"] if row["ontomaId"] else return,
                "cohortPhenotypes" : row["cohortPhenotypes"],
                "targetFromSourceId" : row["Symbol"],
                "allelicRequirements" : row["Mode of inheritance"],
                "studyId" : row["Panel Id"],
                "studyOverview" : row["Panel Name"],
                "literature" : row["publications"]
            }
            return evidence
        except Exception as e:
            logging.error(f'Evidence generation failed for row: {row.name}')
            raise
    
    '''
    @staticmethod
    def removeRedundancies(evidences):
        # Parsing data from the evidences object
        parsed_data = []

        json_data = evidences.apply(lambda x: json.loads(x))
        for evidence in json_data:
            parsed_data.append({
                'target': evidence['target']['id'],
                'disease': evidence['disease']['id'],
                'panel_id': evidence['unique_association_fields']['study_id'],
                'phenotype': evidence['disease']['source_name'],
                'json_data': evidence
            })
        panelapp_df = pd.DataFrame(parsed_data)  
        
        # Grouping to make the evidence unique: by target, disease and panel id
        updated_evidences = []
        for (target, disease, panel_id), group in panelapp_df.groupby(['target','disease','panel_id']):
            # Extracting evidence data:
            data = group["json_data"].tolist()[0]
            # Update evidence:
            source_phenotypes = group.phenotype.unique().tolist()
            data['disease']['source_name'] = choice(source_phenotypes)
            # Save evidence:
            updated_evidences.append(data)

        return updated_evidences
    '''

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates Genomics England PanelApp sourced evidences.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input .tsv file with the table containing association details.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the json output file containing the evidence strings.")
    parser.add_argument("-d", "--mappingsDict", required=False, type=str, help="Path of the dictionary containing the mapped phenotypes.")
    parser.add_argument("-m", "--mappingStep", required=False, type=bool, default=True, help="State whether to run the disease to EFO term mapping step or not.")

    # Parsing parameters
    args = parser.parse_args()

    dataframe = args.inputFile
    outputFile = args.outputFile
    mappingStep = args.mappingStep
    mappingsDict = args.mappingsDict

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidenceBuilder = PanelAppEvidenceGenerator(dataframe, mappingStep, mappingsDict)

    # Import dictionary if present
    if mappingsDict:
        with open(mappingsDict, "r") as f:
                phenotypesMappings = json.load(f)
    else:
        phenotypesMappings = {}
    
    # Writing evidence strings into a json file
    evidences = evidenceBuilder.writeEvidenceFromSource()

    with open(outputFile, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n'))
    
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()