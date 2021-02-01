from common.HGNCParser import GeneParser
import logging
import requests
import argparse
import json
import numpy as np
import gzip
from pyspark import SparkContext, SparkFiles
from pyspark.sql import SparkSession, Row
from pyspark.sql.functions import *
from pyspark.sql.types import *

class phewasEvidenceGenerator():
    def __init__(self):
        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('phewas') \
                .getOrCreate()

        # Initialize mapping variables
        self.mappingStep = mappingStep
        self.mappingsFile = "https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv"

        # Data extracted from Varient Index
        self.consequencesFile = "https://storage.googleapis.com/otar000-evidence_input/PheWAS/data_files/phewas_w_consequences.csv"

        # Adding remote files to Spark Context
        self.spark.sparkContext.addFile(self.consequencesFile)
    
        # Initialize input files
        self.inputFile = inputFile
        
        # Initialize gene parser
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.udfGeneParser = udf(
            lambda X: gene_parser.genes.get(X.strip("*"), np.nan),
            StringType()
        )

        self.dataframe = None
        self.enrichedDataframe = None

    def generateEvidenceFromSource(self, inputFile, diseaseMapping, skipMapping):
        '''
        Processing of the dataframe to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''
        # Read input file
        self.dataframe = self.spark.read.csv(self.inputFile, header=True)

        # Filter out null genes & p-value > 0.05
        self.dataframe = self.dataframe \
                        .filter(col("gene").isNotNull()) \
                        .filter(col("p") < 0.05)

        # Mapping step
        if not skipMapping:
            self.spark.sparkContext.addFile(diseaseMapping)
            phewasMapping = self.spark.read.csv(SparkFiles.get("phewascat.mappings.tsv"), sep=r'\t', header=True)
            self.dataframe = self.dataframe.join(
                phewasMapping,
                on=["Phewas_string"],
                how="inner"
            )
        
        self.dataframe = self.dataframe.withColumn(
            "gene",
            self.udfGeneParser(col("gene"))
        )

        # Get functional consequence per variant from OT Genetics Portal
        cols = ["phewas_string", "phewas_code", "EFO_id", "odds_ratio", "p", "cases", "gene", "consequence_link", "variantId", "snp"]
        self.enrichedDataframe = self.enrichVariantData() \
                                        .dropDuplicates(cols)

        # Build evidence strings per row
        evidences = self.enrichedDataframe.rdd \
            .map(phewasEvidenceGenerator.parseEvidenceString) \
            .collect() # list of dictionaries
        
        return evidences

    def enrichVariantData(self):
        phewasWithConsequences = self.spark.read.csv(SparkFiles.get("phewas_w_consequences.csv"), header=True)
        phewasWithConsequences = phewasWithConsequences \
                                        .withColumnRenamed("rsid", "snp") \
                                        .withColumnRenamed("gene_id", "gene")

        # Enriching dataframe with consequences --> more records due to 1:many associations
        self.dataframe = self.dataframe.join(
            phewasWithConsequences,
            on=["gene", "snp"],
            how="inner"
        )

        # Building variantId: "chrom_pos_ref_alt" of the respective rsId
        # If one rsId has several variants, variantId = none
        one2manyVariants_df = phewasWithConsequences \
                                    .groupBy("snp") \
                                    .agg(count("snp")) \
                                    .filter(col("count(snp)") > 1)
        one2manyVariants = list(one2manyVariants_df.toPandas()["snp"])
        self.enrichedDataframe = self.dataframe.rdd.map(lambda X: phewasEvidenceGenerator.writeVariantId(X, one2manyVariants)).toDF()
        return self.enrichedDataframe

    @staticmethod
    def writeVariantId(row, one2manyVariants):
        rd = row.asDict()
        if row["snp"] not in one2manyVariants:
            rd["variantId"] = "{}_{}_{}_{}".format(row["chrom"], row["pos"], row["ref"], row["alt"]) # TO-DO: fix pos casting
        else:
            rd["variantId"] = None
        new_row = Row(**rd)
        return new_row

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "phewas_catalog",
                "datatypeId" : "genetic_association",
                "diseaseFromSource" : row["phewas_string"],
                "diseaseFromSourceId" : row["phewas_code"],
                "diseaseFromSourceMappedId" : row["EFO_id"].split("/")[-1],
                "oddsRatio" : row["odds_ratio"],
                "resourceScore" : row["p"],
                "studyCases" : row["cases"],
                "targetFromSourceId" : row["gene"].strip("*"),
                "variantFunctionalConsequenceId" : row["consequence_link"].split("/")[-1] if row["consequence_link"] else "SO_0001060",
                "variantId" : row["variantId"],
                "variantRsId" : row["snp"]
            }
            return evidence
        except Exception as e:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences from the PheWAS Catalog data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input .csv file with the table containing association details.")
    parser.add_argument("-c", "--consequencesFile", required=True, type=str, help="Input look-up table containing the variation consequences coming from the Variant Index.")
    parser.add_argument("-d", "--diseaseMapping", required=False, type=str, help="Input look-up table containing the phenotype mappings to an EFO ID.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the compressed json.gz output file containing the evidence strings.")
    parser.add_argument("-s", "--skipMapping", required=False, action="store_true", help="State whether to skip the disease to EFO mapping step.")
    parser.add_argument("-l", "--logFile", help="Destination of the logs generated by this script.", type=str, required=False)

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    consequencesFile = args.consequencesFile
    diseaseMapping = args.diseaseMapping
    outputFile = args.outputFile
    skipMapping = args.skipMapping

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)
    
    # Logging parameters
    logging.info(f"PheWAS input table: {inputFile}")
    logging.info(f"Phewas phenotype to EFO ID table: {diseaseMapping}")
    logging.info(f"Phewas enriched with consequences input file: {consequencesFile}")
    logging.info(f"Output file: {outputFile}")

    # Initialize evidence builder object
    evidenceBuilder = phewasEvidenceGenerator(dataframe, mappingStep)

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.generateEvidenceFromSource()
        
    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()
    
