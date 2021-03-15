from ontoma import OnToma
import logging
from sys import stderr
import requests
import argparse
import re
import gzip
import json
import multiprocessing as mp # from multiprocessing import Pool
import numpy as np
import pandas as pd
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *
import time

class PanelAppEvidenceGenerator():
    def __init__(self, phenotypesMappings):
        # Create OnToma object
        self.otmap = OnToma()

        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('evidence_builder') \
                .getOrCreate()

        self.dataframe = None
        
        # Initialize mapping variables
        self.diseaseMappings = phenotypesMappings
        self.codesMappings = None

    def writeEvidenceFromSource(self, inputFile, skipMapping):
        '''
        Processing of the dataframe to build all the evidences from its data

        Args:
            dataframe (pyspark.DataFrame): Initial .tsv file
        Returns:
            evidences (array): Object with all the generated evidences strings
        '''

        # Reading and filtering input file
        self.dataframe = (self.spark
                            .read.csv(inputFile, sep=r'\t', header=True)
                            .filter(
                                ((col("List") == "green" ) | (col("List") == "amber")) &
                                (col("Panel Version") > 1) &
                                (col("Panel Status") == "PUBLIC")
                            )
        )

        logging.info("Fetching publications from the API...")
        pdf = PanelAppEvidenceGenerator.buildPublications(self.dataframe.toPandas()) # TODO: write in pyspark 
        pdf.dropna(axis=1, how='all', inplace=True)
        self.dataframe = self.spark.createDataFrame(pdf)
        logging.info("Publications loaded.")

        # Cleaning the phenotype related data of the dataframe
        self.dataframe = PanelAppEvidenceGenerator.cleanDataframe(self.dataframe)

        # Map the diseases to an EFO term if necessary
        if not skipMapping:
            self.dataframe = self.diseaseMappingStep()
        else:
            logging.info("Disease mapping has been skipped.")
            self.dataframe = self.dataframe.withColumn(
                "ontomaUrl",
                lit(None)
            )

        logging.info("Generating evidence:")
        evidences = (self.dataframe
                        # Removing redundant evidence after the explosion of phenotypes
                        .dropDuplicates(["Panel Id", "Symbol", "ontomaUrl"])
                        # Build evidence strings per row
                        .rdd.map(
                            PanelAppEvidenceGenerator.parseEvidenceString)
                        .collect() # list of dictionaries
        )
        
        for evidence in evidences:
            # Delete empty keys
            if skipMapping:
                del evidence["diseaseFromSourceMappedId"]
            if not evidence["diseaseFromSourceId"]:
                del evidence["diseaseFromSourceId"]
            if len(evidence["literature"]) == 0:
                del evidence["literature"]

        return evidences

    @staticmethod
    def buildPublications(pdf):
        '''
        Populates a dataframe with the publications fetched from the PanelApp API and cleans them to match PubMed IDs.

        Args:
            dataframe (pandas.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pandas.DataFrame): DataFrame with an "publications" column added
        '''
        populated_groups = []

        for (PanelId), group in pdf.groupby("Panel Id"):
            request = PanelAppEvidenceGenerator.publicationsFromPanel(PanelId)
            group["publications"] = group.apply(lambda X: PanelAppEvidenceGenerator.publicationFromSymbol(X.Symbol, request), axis=1)
            populated_groups.append(group)
        
        pdf = pd.concat(populated_groups, ignore_index=True, sort=False)

        cleaned_publication = []
        for row in pdf["publications"].to_list():
            try:
                tmp = [re.match(r"(\d{8})", e)[0] for e in row]
                cleaned_publication.append(list(set(tmp))) # Removing duplicated publications
            except Exception as e:
                cleaned_publication.append([])
                continue

        pdf["publications"] = cleaned_publication

        return pdf

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
            logging.error("Query of the PanelApp API has failed.")
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
    
        dataframe = (dataframe
                        .withColumn(
                            # NaNs and "No OMIM phenotype" in "Phenotypes" column --> Assignment of Panel Name
                            "Phenotypes",
                            when(col("Phenotypes") == "No OMIM phenotype",
                                col("Panel Name"))
                            .when(col("Phenotypes").isNull(),
                                col("Panel Name"))
                            .otherwise(col("Phenotypes"))
                        )
                        # cohortPhenotypes --> array of the original string separated by phenotypes 
                        .withColumn(
                            "cohortPhenotypes",
                            array_distinct(split(
                                col("Phenotypes"),
                                ";"
                            ))   
                        )
                        # phenotype --> explosion of cohortPhenotypes
                        .withColumn(
                            "phenotype",
                            explode(col("cohortPhenotypes"))
                        )
                        # omimCode --> OMIM code present in "phenotype"
                        .withColumn("omimCode", regexp_extract(col("phenotype"), "(\d{6})", 1))
                        # removal of the OMIM code in "phenotype"
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "(\d{6})", ""))
                        # deleting special characters in "phenotype"
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "[^0-9a-zA-Z -]", ""))
                        .withColumn("phenotype", trim(col("phenotype")))
        )

        return dataframe
    
    def diseaseMappingStep(self):
        '''
        Runs the disease mapping step from phenotypes and OMIM codes to an EFO ID - 3 steps:
            1. Querying OnToma with all distinct codes and phenotypes
            2. Xref between the results of every phenotype/OMIM code pair for a better coverage
            3. Mappings are built into the dataframe

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pyspark.DataFrame): DataFrame with the mapping results filtered by only matches
        '''

        omimCodesDistinct = list(self.dataframe.select("omimCode").distinct().toPandas()["omimCode"])
        phenotypesDistinct = list(self.dataframe.select("phenotype").distinct().toPandas()["phenotype"])

        if self.diseaseMappings is None:
            # Checks whether the dictionary is not provided as a parameter
            try:
                pool = mp.Pool(mp.cpu_count())
                self.diseaseMappings = pool.map(self.diseaseToEfo, phenotypesDistinct) # list of dictionaries
                self.diseaseMappings = {k:v for dct in self.diseaseMappings for (k,v) in self.diseaseMappings.items()}
                pool.close()
            except Exception as e:
                logging.error(f"Processing the mappings without multithreading. An error occurred during the task: \n{e}.")
                self.diseaseMappings = self.diseaseToEfo(*phenotypesDistinct)
        else:
            logging.info(f"Disease mappings have been imported.")

        self.codesMappings = self.diseaseToEfo(*omimCodesDistinct, dictExport="codesToEfo_results.json")
        for i, (phenotype, code) in self.dataframe.select("omimCode", "phenotype").toPandas().drop_duplicates().iterrows():
            self.phenotypeCodePairCheck(phenotype, code)
        logging.info("Disease mappings have been checked.")

        # Add new columns: ontomaResult, ontomaUrl, ontomaLabel
        phenotypesMappings = self.diseaseMappings
        udfBuildMapping = udf(
            lambda X: PanelAppEvidenceGenerator.buildMapping(X, phenotypesMappings),
            ArrayType(StringType()))
        self.dataframe = (self.dataframe
            .withColumn("ontomaResult", udfBuildMapping(col("phenotype"))[0])
            .withColumn("ontomaUrl",
                element_at(
                    split(
                        udfBuildMapping(col("phenotype"))[1],
                        "/"),
                    -1
                )
            )
            .withColumn("ontomaLabel", udfBuildMapping(col("phenotype"))[1])
        )

        return self.dataframe.filter(col("ontomaResult") == "match")

    def diseaseToEfo(self, *iterable, dictExport="diseaseToEfo_results.json"):
        '''
        Queries the OnToma utility to map a phenotype to an EFO ID.

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
        Among the Fuzzy results of a phenotype query, it checks if the phenotype and the respective code points
        to the same EFO term and changes the dictionary with the mappings

        Args:
            phenotype (str): Phenotype of each phenotype-OMIM code pair
            omimCode (str): Code of each phenotype-OMIM code pair
        '''
        try:
            if phenotype == "":
                return
            phenotypeResult = self.diseaseMappings[phenotype]
            if phenotypeResult == None:
                return

            if phenotypeResult["quality"] == "fuzzy":
                codeResult = self.codesMappings[code]
                if codeResult == None:
                    return

                if codeResult["term"] == phenotypeResult["term"]:
                    # If both EFO terms coincide, the phenotype in mappingsDict becomes a match
                    self.diseaseMappings[phenotype]["quality"] = "match"
                    self.diseaseMappings[phenotype]["action"] = "checked"
        except:
            logging.error(f'No OMIM code for phenotype: {phenotype}')
    
    def buildMapping(phenotype, phenotypesMappings):
        if phenotype in phenotypesMappings.keys():
            ontomaResult = phenotypesMappings[phenotype]["quality"]
            ontomaUrl = phenotypesMappings[phenotype]["term"]
            ontomaLabel = phenotypesMappings[phenotype]["label"]

            return ontomaResult, ontomaUrl, ontomaLabel # TO-DO: implement this https://stackoverflow.com/questions/42980704/pyspark-create-new-column-with-mapping-from-a-dict

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "genomics_england",
                "datatypeId" : "genetic_literature",
                "confidence" : row["List"],
                "diseaseFromSource" : row["phenotype"],
                "diseaseFromSourceId" : row["omimCode"],
                "diseaseFromSourceMappedId" : row["ontomaUrl"],
                "cohortPhenotypes" : row["cohortPhenotypes"],
                "targetFromSourceId" : row["Symbol"],
                "allelicRequirements" : [row["Mode of inheritance"]],
                "studyId" : row["Panel Id"],
                "studyOverview" : row["Panel Name"],
                "literature" : row["publications"]
            }
            return evidence
        except Exception as e:
            logging.error(f'Evidence generation failed for row: {"row"}')
            raise

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates Genomics England PanelApp sourced evidences.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input .tsv file with the table containing association details.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the compressed json.gz output file containing the evidence strings.")
    parser.add_argument("-d", "--mappingsDict", required=False, type=str, help="JSON file containing the mapped phenotypes of an earlier run.")
    parser.add_argument("-s", "--skipMapping", required=False, action="store_true", help="State whether to skip the disease to EFO mapping step.")
    parser.add_argument("-l", "--logFile", help="Destination of the logs generated by this script.", type=str, required=False)

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    outputFile = args.outputFile
    skipMapping = args.skipMapping
    mappingsDict = args.mappingsDict

    # Initialize logging:
    logging.basicConfig(
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(stderr)

    # Logging parameters
    logging.info(f"Panel app input table: {inputFile}")
    logging.info(f"Phenotypes mapped to EFO look-up file: {mappingsDict}")
    logging.info(f"Output file: {outputFile}")

    # Import dictionary if present
    if mappingsDict:
        with open(mappingsDict, "r") as f:
                phenotypesMappings = json.load(f)
    else:
        phenotypesMappings = None

    # Initialize evidence builder object
    evidenceBuilder = PanelAppEvidenceGenerator(phenotypesMappings)

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.writeEvidenceFromSource(inputFile, skipMapping)

    # Exporting the outfile
    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"{len(evidences)} evidence strings have been generated.")

if __name__ == '__main__':
    main()