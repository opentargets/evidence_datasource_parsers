#!/usr/bin/env python

import sys
import argparse
import gzip
import logging
import json
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

class progenyEvidenceGenerator():
    def __init__(self):
        # Create spark session     
        self.spark = (SparkSession.builder
                .appName('progeny')
                .getOrCreate())
    
        # Initialize source table
        self.dataframe = None

    def generateEvidenceFromSource(self, inputFile, diseaseMapping, pathwayMapping, skipMapping):
        '''
        Processing of the input file to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''
        # Read input file
        self.dataframe = (self.spark.read
                                .option("header", "true")
                                .option("delimiter", "\t")
                                .option("inferSchema", "true")
                                .csv(inputFile))

        # Disease mapping step
        if not skipMapping:
            try:
                self.dataframe = self.cancer2EFO(diseaseMapping)
                logging.info("Disease mappings have been imported.")
            except Exception as e:
                logging.error(f"An error occurred while importing disease mappings: \n{e}.")
        else:
            logging.info("Disease mapping has been skipped.")
            self.dataframe = self.dataframe.withColumn(
                "EFO_id",
                lit(None)
            )
        
        self.dataframe = self.pathway2Reactome(pathwayMapping)
        logging.info("Pathway to reaction ID mappings have been imported.")

        # Build evidence strings per row
        logging.info("Generating evidence:")
        evidences = (self.dataframe.rdd
            .map(progenyEvidenceGenerator.parseEvidenceString)
            .collect() # list of dictionaries
            )

        if skipMapping:
            # Delete empty keys if mapping is skipped
            for evidence in evidences:
                del evidence["diseaseFromSourceMappedId"]

        return evidences
    
    def cancer2EFO(self, diseaseMapping):
        diseaseMappingsFile = (self.spark
                        .read.csv("resources/cancer2EFO_mappings.tsv", sep=r'\t', header=True)
                        .select("Cancer_type_acronym", "EFO_id")
                        .withColumnRenamed("Cancer_type_acronym", "Cancer_type"))

        self.dataframe = self.dataframe.join(
            diseaseMappingsFile,
            on="Cancer_type",
            how="inner"
        )

        return self.dataframe
    
    def pathway2Reactome(self, pathwayMapping):
        pathwayMappingsFile = (self.spark
                        .read.csv("resources/pathway2Reactome_mappings.tsv", sep=r'\t', header=True)
                        .withColumnRenamed("pathway", "Pathway"))
        
        self.dataframe = (self.dataframe
                .join(
                    pathwayMappingsFile,
                    on="Pathway",
                    how="inner"
                )
                .withColumn(
                    "target",
                    split(col("target"), ", ") 
                )
                .withColumn(
                    "target",
                    explode("target")
                ))

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "progeny",
                "datatypeId" : "affected_pathway",
                "diseaseFromSourceMappedId" : row["EFO_id"],
                "resourceScore" : row["P.Value"],
                "targetFromSourceId" : row["target"],
                "diseaseFromSource" : row["Cancer_type"],
                "pathways" : [
                    {
                    "id" : row["reactomeId"],
                    "name" : row["description"]
                    }
                ]
            }
            return evidence
        except Exception:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences for the PROGENy data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input source .tsv file.")
    parser.add_argument("-d", "--diseaseMapping", required=False, type=str, help="Input look-up table containing the cancer type mappings to an EFO ID.")
    parser.add_argument("-p", "--pathwayMapping", required=True, type=str, help="Input look-up table containing the pathway mappings to a respective target and ID in Reactome.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Gzipped JSON file containing the evidence strings.")
    parser.add_argument("-s", "--skipMapping", required=False, action="store_true", help="State whether to skip the disease to EFO mapping step.")
    parser.add_argument("-l", "--logFile", help="Destination of the logs generated by this script.", type=str, required=False)

    # Parsing parameters
    args = parser.parse_args()
    inputFile = args.inputFile
    diseaseMapping = args.diseaseMapping
    pathwayMapping = args.pathwayMapping
    outputFile = args.outputFile
    skipMapping = args.skipMapping

    # Initialize logging:
    logging.basicConfig(
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    # Logging parameters
    logging.info(f"PROGENy input table: {inputFile}")
    logging.info(f"Cancer type to EFO ID table: {diseaseMapping}")
    logging.info(f"Pathway to reaction ID table: {pathwayMapping}")
    logging.info(f"Output file: {outputFile}")

    # Initialize evidence builder object
    evidenceBuilder = progenyEvidenceGenerator()

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.generateEvidenceFromSource(inputFile, diseaseMapping, pathwayMapping, skipMapping)

    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"{len(evidences)} evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()