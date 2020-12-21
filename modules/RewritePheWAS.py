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
        self.mappingsFile = "phewascat.mappings.tsv" # TODO : Remove this patch
        phewasMapping = spark.read.csv(mappingFile, sep=r'\t', header=True)

        # Filter out null genes & p-value > 0.05
        dataframe = dataframe \
                        .filter(col("gene").isNull() == False) \
                        .filter(col("p") < 0.05)

        # Mapping step
        if self.mappingStep:
            dataframe = dataframe.join(
                phewasMapping,
                on=["Phewas_string"],
                how="inner"
            )
        # Build evidence strings per row
        evidences = dataframe.rdd.map(
            phewasEvidenceGenerator.parseEvidenceString)
            .collect() # list of dictionaries
        
        pass
    
    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                datasourceId : "phewas_catalog",
                datatypeId : "genetic_association",
                diseaseFromSource : row["phewas_string"],
                diseaseFromSourceId : row["phewas_code"],
                diseaseFromSourceMappedId : row["EFO_id"].split("/")[-1],
                oddsRatio : row["odds_ratio"],
                resourceScore : row["p"],
                studyCases : row["cases"],
                targetFromSourceId : row["gene"].strip("*"),
                variantFunctionalConsequenceId : "", # TODO : Merge this data coming from OTG
                variantId : "", # TODO 
                variantRsId : row["snp"]
            }
            return evidence
        except Exception as e:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences from the PheWAS Catalog data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input .csv file with the table containing association details.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the json output file containing the evidence strings.")
    parser.add_argument("-s", "--schemaVersion", required=True, type=str, help="JSON schema version to use, e.g. 1.6.9. It must be branch or a tag available in https://github.com/opentargets/json_schema.")
    parser.add_argument("-m", "--mappingStep", required=True, type=bool, default=True, help="State whether to run the disease to EFO term mapping step or not.")

    # Parsing parameters
    args = parser.parse_args()

    dataframe = args.inputFile
    outputFile = args.outputFile
    schemaVersion = args.schemaVersion
    mappingStep = args.mappingStep

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',E
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidenceBuilder = phewasEvidenceGenerator(schemaVersion, mappingStep)

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.writeEvidenceFromSource(dataframe)

    with open(outputFile, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n'))
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()
    pass
