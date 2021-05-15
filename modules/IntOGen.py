#!/usr/bin/env python

import argparse
import sys
import gzip
import logging
import json

from pyspark.sql import SparkSession
import pyspark.sql.functions as F

class intogenEvidenceGenerator():
    def __init__(self):
        
        # Create spark session
        self.spark = (
            SparkSession.builder
            .appName('intOGen')
            .getOrCreate()
        )
        
        logging.info(f"Spark version: {self.spark.version}")

        # Initialize source tables
        self.dataframe = None

    def generateEvidenceFromSource(self, inputGenes, inputCohorts, diseaseMapping, skipMapping):
        '''
        Processing of the input file to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''

        genes = (
            self.spark
            .read.csv(inputGenes, sep=r'\t', header=True, inferSchema=True)
            .select(
                'SYMBOL',
                'COHORT',
                F.col("CANCER_TYPE").alias("Cancer_type_acronym"),
                F.col("SAMPLES").alias("numberMutatedSamples"),
                "METHODS",
                "ROLE",
                "QVALUE_COMBINATION"
            )
            .withColumn("METHODS", F.split(F.col("METHODS"), ","))
          
            # Mutation role mapping to a SO code
            .withColumn(
                "functionalConsequenceId",
                # gain_of_function_variant
                F.when(F.col("ROLE") == "Act", "SO_0002053")
                # loss_of_function_variant
                .when(F.col("ROLE") == "LoF", "SO_0002054")
                .otherwise(None)
            )
        )

        cohorts = (
            self.spark
            .read.csv(inputCohorts, sep=r'\t', header=True, inferSchema=True)
            .select(
                "COHORT",
                "CANCER_TYPE_NAME",
                "WEB_SHORT_COHORT_NAME",
                "WEB_LONG_COHORT_NAME",
                F.col("SAMPLES").alias("numberSamplesTested")
            )
        )

        # Joining genes and cohorts data
        self.dataframe = genes.join(
            cohorts,
            on="COHORT",
            how="inner"
        )

        # Mapping step
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
                F.lit(None)
            )

        # Build evidence strings per row
        logging.info("Generating evidence:")
        evidences = (
            self.dataframe.rdd
            .map(intogenEvidenceGenerator.parseEvidenceString)
            .collect()
        )  # list of dictionaries

        for evidence in evidences:
            if skipMapping:
                # Delete empty keys if mapping is skipped
                del evidence["diseaseFromSourceMappedId"]
            if not evidence["mutatedSamples"][0]["functionalConsequenceId"]:
                # Dele empty key if functionalConsequenceId is unknown
                del evidence["mutatedSamples"][0]["functionalConsequenceId"]

        return evidences

    def cancer2EFO(self, diseaseMapping):

        diseaseMappingsFile = (
            self.spark
            .read.csv(diseaseMapping, sep=r'\t', header=True)
            .select("Cancer_type_acronym", "EFO_id")
            .withColumn("EFO_id", F.trim(F.col("EFO_id")))
        )

        self.dataframe = self.dataframe.join(
            diseaseMappingsFile,
            on="Cancer_type_acronym",
            how="left"
        )

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId": "intogen",
                "cohortDescription": row["WEB_LONG_COHORT_NAME"],
                "cohortId": row["COHORT"],
                "cohortShortName": row["WEB_SHORT_COHORT_NAME"],
                "datatypeId": "somatic_mutation",
                "diseaseFromSource": row["CANCER_TYPE_NAME"],
                "resourceScore": row["QVALUE_COMBINATION"],
                "significantDriverMethods": row["METHODS"],
                "targetFromSourceId": row["SYMBOL"],
                "mutatedSamples": [{
                    "functionalConsequenceId": row["functionalConsequenceId"],
                    "numberMutatedSamples": row["numberMutatedSamples"],
                    "numberSamplesTested": row["numberSamplesTested"]
                }]
            }

            # The diseaseFromSourceMappedId is skipped if EFO_id is missing:
            if row['EFO_id']:
                evidence['diseaseFromSourceMappedId'] = row['EFO_id']

            return evidence
        except Exception as e:
            raise


def main(inputGenes, inputCohorts, diseaseMapping, outputFile, skipMapping):

    # Logging parameters
    logging.info(f"intOGen driver genes table: {inputGenes}")
    logging.info(f"intOGen cohorts table: {inputCohorts}")
    logging.info(f"Cancer type to EFO ID table: {diseaseMapping}")
    logging.info(f"Output file: {outputFile}")
    logging.info(f"Skipping disease mapping: {skipMapping}")

    # Initialize evidence builder object
    evidenceBuilder = intogenEvidenceGenerator()

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.generateEvidenceFromSource(inputGenes, inputCohorts, diseaseMapping, skipMapping)

    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')

    logging.info(f"{len(evidences)} evidence strings saved into {outputFile}. Exiting.")


if __name__ == '__main__':

    # Initiating parser
    parser = argparse.ArgumentParser(description="This script generates evidences for the intOGen data source.")

    parser.add_argument("-g", "--inputGenes", required=True, type=str, help="Input source .tsv file listing the driver genes across the analyzed cohorts.")
    parser.add_argument("-c", "--inputCohorts", required=True, type=str, help="Input source .tsv file with information about the analyzed cohorts.")
    parser.add_argument("-d", "--diseaseMapping", required=False, type=str, help="Input look-up table containing the cancer type mappings to an EFO ID.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Gzipped JSON file containing the evidence strings.")
    parser.add_argument("-s", "--skipMapping", required=False, action="store_true", help="State whether to skip the disease to EFO mapping step.")
    parser.add_argument("-l", "--logFile", help="Destination of the logs generated by this script.", type=str, required=False)

    # Parsing parameters
    args = parser.parse_args()
    inputGenes = args.inputGenes
    inputCohorts = args.inputCohorts
    diseaseMapping = args.diseaseMapping
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

    main(inputGenes, inputCohorts, diseaseMapping, outputFile, skipMapping)
