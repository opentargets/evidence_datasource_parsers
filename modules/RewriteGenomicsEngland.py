from common.HGNCParser import GeneParser
from ontoma import OnToma
import logging
import requests
import argparse
import re
from pyspark import SparkContext
from pyspark.sql import SparkSession
import pyspark.sql.functions import col, coalesce, when, udf, explode, regexp_extract, regexp_replace
from pyspark.sql.types import *

class PanelEvidenceGenerator():
    def __init__(self, schema_version):
        # Build JSON schema url from version
        self.schema_version = schema_version
        schema_url = f"https://raw.githubusercontent.com/opentargets/json_schema/{self.schema_version}/draft4_schemas/opentargets.json" #TODO Update the url 
        logging.info(f"Loading JSON Schema from {schema_url}")

        # Create OnToma object
        self.otmap = OnToma()

        # Parse Gene Symbol to EnsEMBL ID
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('evidence_builder') \
                .getOrCreate()

    def writeEvidenceFromSource(self, dataframe):
        '''
        Processing of the dataframe to build all the evidences from its data

        Args:
            dataframe (pd.DataFrame): Initial .tsv file
            mappings_dict (dict): All mapping results for every phenotype
        Returns:
            evidences (array): Object with all the generated evidences strings
        '''

        # Read input file
        dataframe = self.spark.read.csv(dataframe, sep=r'\t', header=True)

        # Filtering with green and ambar evidences
        dataframe = dataframe.filter(
                                ((col("List") == "green" ) | (col("List") == "amber"))
                                                            &
                                (col("Panel Version") > 1)
                                                            &
                                (col("Panel Status") == "PUBLIC")                
        )

        # TODO: Feed the dataframe with publications
        logging.info("Fetching publications from the API...")
        #dataframe = self.build_publications(dataframe)
        logging.info("Publications loaded.")

        # Splitting and cleaning the dataframe according to the phenotype string
        dataframe = PanelAppEvidenceGenerator.clean_dataframe(dataframe)
        pass

    def cleanDataframe(dataframe):
        '''
        Args:
            dataframe (pd.DataFrame): Initial .tsv data converted to a Pandas dataframe
        Returns:
            dataframe (pd.DataFrame): Original dataframe cleaned
        '''
        # TODO: Wrap different steps into one function to iterate less times
        
        # NaNs and "No OMIM phenotype" in "Phenotypes" column --> Assignment of Pannel Name
        dataframe = dataframe.withColumn("cohortPhenotypes", when(col("Phenotypes") == "No OMIM phenotype", col("Panel Name"))
                                 .when(col("Phenotypes").isNull(), col("Panel Name"))
                                 .otherwise(col("Phenotypes")))
        
        # Handling multiple phenotypes column: 
        #   cohortPhenotypes --> original string 
        #   phenotype --> explosion of cohortPhenotypes
        splitLambda = udf(lambda X: X.split(";"), ArrayType(StringType()))
        dataframe = dataframe \
                            .withColumn("phenotype", splitLambda(col("cohortPhenotypes"))) \
                            .withColumn('phenotype', explode('phenotype'))

        # Extracting and cleaning the OMIM codes: 
        #   removal of the OMIM codes in the Phenotypes column and the inclusion in omimCodes
        #   deleting special characters
        stripLambda = udf(lambda X: X.strip(),StringType())
        dataframe = dataframe \
                        .withColumn("omimCode", regexp_extract(col("phenotype"), "(\d{6})", 1)) \
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "(\d{6})", "")) \
                        .withColumn("phenotype", regexp_replace(col("phenotype"), "[^0-9a-zA-Z *]", "")) \
                        .withColumn("phenotype", stripLambda(col("phenotype")))



        pass
