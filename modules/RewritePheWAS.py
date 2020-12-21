import logging
import requests
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

class phewasEvidenceGenerator():
    def __init__(self, schemaVersion, mappingStep):
        # Build JSON schema url from version
        self.schemaVersion = schemaVersion
        schema_url = f"https://raw.githubusercontent.com/opentargets/json_schema/{self.schemaVersion}/draft4_schemas/opentargets.json" #TODO Update the url 
        logging.info(f"Loading JSON Schema from {schema_url}")

        # Initialize mapping variables
        self.mappingStep = mappingStep
        self.mappingsFile = "https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv" 

        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('evidence_builder') \
                .getOrCreate()

    def writeEvidenceFromSource(self, dataframe):
        '''
        Processing of the dataframe to build all the evidences from its data
        Args:
            dataframe (.csv): Source file containing all data from PheWAS
        Returns:
            evidences (array): Object with all the generated evidences strings
        '''
        
        # Read input file + manually mapped terms 
        dataframe = self.spark.read.csv("")
        mappingFile = "phewascat.mappings.tsv" # TODO : Should be "https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv" but it is failing
        phewasMapping = spark.read.csv(mappingFile, sep=r'\t', header=True)

        # Filter out null genes & p-value > 0.05
        dataframe = dataframe \
                        .filter(col("gene").isNull() == False) \
                        .filter((col("p") < 0.05)

        # Join mappings data
         dataframe = dataframe.join(
            phewasMapping,
            on=["Phewas_string"],
            how="inner"
         )

        pass
        
