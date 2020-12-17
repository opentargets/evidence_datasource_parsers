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

class PanelAppEvidenceGenerator():
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

    def writeEvidenceFromSource(self, dataframe, phenotypesMappings):
        '''
        Processing of the dataframe to build all the evidences from its data

        Args:
            dataframe (pyspark.DataFrame): Initial .tsv file
            phenotypesMappings (dict): All mapping results for every phenotype
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
        dataframe = PanelAppEvidenceGenerator.cleanDataframe(dataframe)

        # Mapping the phenotypes and OMIM codes to an EFO code - 2 steps:
        #   Querying OnToma with all distinch codes and phenotypes
        #   Xref between the results of every phenotype/OMIM code pair for a better coverage
        omimCodesDistinct = dataframe.select("omimCode").distinct().rdd.flatMap(lambda x: x).collect()
        phenotypesDistinct = dataframe.select("phenotype").distinct().rdd.flatMap(lambda x: x).collect()
        omimCodes = dataframe.select("omimCode").rdd.flatMap(lambda x: x).collect()
        phenotypes = dataframe.select("phenotype").rdd.flatMap(lambda x: x).collect()
        phenotypeCodePairs = list(zip(phenotypes, omimCodes))

        if len(phenotypesMappings) == 0:
            # Checks whether the dictionary is not provided as a parameter 
            phenotypesMappings = self.diseaseToEfo(phenotypesDistinct)
            logging.info("Disease mappings completed.")
        else:
            logging.info("Disease mappings imported.")
        codesMappings = self.diseaseToEfo(omimCodesDistinct) # TODO: add posibility to provide dict

        for pheno, code in phenotypeCodePairs:
            PanelAppEvidenceGenerator.phenotypeCodePairCheck(pheno, code, phenotypesMappings, codesMappings)

        # TODO: Add new columns: OnToma Result, OnToma Term, OnToma Label
        #dataframe = PanelAppEvidenceGenerator.buildMappings(dataframe, phenotypesMappings)

        # Build evidence strings per row
        dataframe = dataframe.filter(col("ontomaResult") == "match")
        evidences = dataframe.apply(self.get_evidence_string, axis=1)
        logging.info(f"{len(evidences)} evidence strings have been generated.")

        pass
    
    @staticmethod
    def cleanDataframe(dataframe):
        '''
        Args:
            dataframe (pyspark.DataFrame): Initial .tsv data converted to a Pandas dataframe
        Returns:
            dataframe (pyspark.DataFrame): Original dataframe cleaned
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

        return dataframe


    def diseaseToEfo(self, iterable, dictExport="diseaseToEfo_results.json"):
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

    @staticmethod
    def phenotypeCodePairCheck(phenotype, omimCode, phenotypesMappings, codesMappings):
        '''
        Among the Fuzzy results of a phenotype query, it checks if the phenotype and the respective code points to the same EFO term

        Args:
            phenotype (str): Phenotype of each phenotype-OMIM code pair
            omimCode (str): Code of each phenotype-OMIM code pair
            phenotypesMappings (dict): All mapping results for every phenotype
            codesMappings (dict): All mapping results for every OMIM code
        '''
        try:
            if phenotype == "":
                return
            phenotypeResult = phenotypesMappings[phenotype]
            if phenotypeResult == None:
                return

            if phenotypeResult["quality"] == "fuzzy":
                codeResult = codesMappings[code]
                if codeResult == None:
                    return

                if codeResult["term"] == phenotypeResult["term"]:
                    # If both EFO terms coincide, the phenotype in mappings_dict becomes a match
                    phenotypesMappings[phenotype]["quality"] = "match"
                    phenotypesMappings[phenotype]["action"] = "checked"
        except Exception as e:
            logging.error(f'No OMIM code for phenotype: {phenotype}')

    @staticmethod
    def buildMappings(dataframe, phenotypesMappings):
        '''
        Populates the dataframe with the mappings resulted from OnToma.

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
            phenotypesMappings (dict): All mapping results for every phenotype
        Return:
            dataframe (pyspark.DataFrame): DataFrame with new columns corresponding to the OnToma result
        '''

        pass

    def parseEvidenceString(self, row):
        # Association level information
        self.targetFromSourceId = row["Symbol"]
        self.diseaseFromSource = row["phenotype"]
        self.diseaseFromSourceId = row["omimCode"]
        self.diseaseId = row["ontomaId"]
        if self.diseaseId is None or self.diseaseId == '':
            logging.warning(f'No EFO id for association row: {row.name}')
            return None
        self.cohortPhenotypes = row["cohortPhenotypes"]
        # self.mapped_disease = row["ontomaLabel"] This field?
        self.allelicRequirements = row["Mode of inheritance"]
        self.confidence = row["List"]
        self.literature = row["publications"]
        try:
            self.targetId = self.genes[self.gene_symbol]
        except:
            self.targetId = row["EnsemblId(GRch37)"]
        self.studyId = row["Panel Id"]
        self.studyOverview = row["Panel Name"]

        try:
            evidence = {
                "datasourceId" : "genomics_england",
                "datatypeId" : "genetic_literature",
                "confidence" : self.confidence,
                "diseaseFromSource" : self.diseaseFromSource,
                "diseaseFromSourceId" : self.diseaseFromSourceId,
                "diseaseId" : self.diseaseId,
                "cohortPhenotypes" : self.cohortPhenotypes,
                "targetFromSourceId" : self.targetFromSourceId,
                "allelicRequirements" : self.allelicRequirements,
                "studyId" : self.studyId,
                "studyOverview" : self.studyOverview,
                "literature" : self.literature,
            }
            return evidence.serialize()
        except Exception as e:
            print(e)
            logging.error(f'Evidence generation failed for row: {row.name}')
            raise