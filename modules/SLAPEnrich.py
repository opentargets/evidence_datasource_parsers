import argparse
import gzip
import logging
import json
from pyspark import SparkContext, SparkFiles
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

# *** Tumor acronym map ***
TUMOR_TYPE_EFO_MAP = {
    'ALL': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000220', 'label': 'acute lymphoblastic leukemia'},
    'BLCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000292', 'label': 'bladder carcinoma'},
    'BRCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000305', 'label': 'breast carcinoma'},
    'CLL': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000095', 'label': 'chronic lymphocytic leukemia'},
    'DLBC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000403', 'label': 'diffuse large B-cell lymphoma'},
    'ESCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0002916', 'label': 'esophageal carcinoma'},
    'GBM': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000519', 'label': 'glioblastoma multiforme'},
    'HNSC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000181', 'label': 'head and neck squamous cell carcinoma'},
    'KIRC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000349', 'label': 'clear cell renal carcinoma'},
    'LAML': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000222', 'label': 'acute myeloid leukemia'},
    'LGG': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0005543', 'label': 'brain glioma'},
    'LIHC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000182', 'label': 'hepatocellular carcinoma'},
    'LUAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000571', 'label': 'lung adenocarcinoma'},
    'LUSC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000708', 'label': 'squamous cell lung carcinoma'},
    'MB': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0002939', 'label': 'medulloblastoma'},
    'MM': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0001378', 'label': 'multiple myeloma'},
    'NB': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000621', 'label': 'neuroblastoma'},
    'OV': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0002917', 'label': 'ovarian serous adenocarcinoma'},
    'PAAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_1000044', 'label': 'pancreatic adenocarcinoma'},
    'PRAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000673', 'label': 'prostate adenocarcinoma'},
    'SCLC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000702', 'label': 'small cell lung carcinoma'},
    'SKCM': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000389', 'label': 'cutaneous melanoma'},
    'STAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000503', 'label': 'stomach adenocarcinoma'},
    'THCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0002892', 'label': 'thyroid carcinoma'},
    'UCEC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_1000233', 'label': 'endometrial endometrioid adenocarcinoma'}
}

# *** Cancer acronyms ***
TUMOR_TYPE_MAP = {
    'ALL': 'acute lymphocytic leukemia',
    'BLCA': 'bladder carcinoma',
    'BRCA': 'breast carcinoma',
    'CLL': 'chronic lymphocytic leukemia',
    'DLBC': 'diffuse large B cell lymphoma',
    'ESCA': 'esophageal carcinoma',
    'GBM': 'glioblastoma multiforme',
    'HNSC': 'head and neck squamous cell carcinoma',
    'KIRC': 'clear cell renal carcinoma',
    'LAML': 'acute myeloid leukemia',
    'LGG': 'lower grade glioma',
    'LIHC': 'hepatocellular carcinoma',
    'LUAD': 'lung adenocarcinoma',
    'LUSC': 'lung squamous cell carcinoma',
    'MB': 'medulloblastoma',
    'MM': 'multiple myeloma',
    'NB': 'neuroblastoma',
    'OV': 'serous ovarian adenocarcinoma',
    'PAAD': 'pancreas adenocarcinoma',
    'PRAD': 'prostate adenocarcinoma',
    'SCLC': 'small cell lung carcinoma',
    'SKCM': 'cutaneous melanoma',
    'STAD': 'stomach adenocarcinoma',
    'THCA': 'thyroid carcinoma',
    'UCEC': 'endometrial endometrioid adenocarcinoma'
}

# *** TO DO: Symbol mapping for genes withouth Ensembl IDs ***
SYMBOL_MAPPING = {
    # TODO These symbols do not have Ensembl ID mappings, need alternative mapping
    '''
    i.e 'C15orf55': 'NUTM1',

    ADRBK1
    ADRBK2
    BRE
    C10orf2
    CASC5
    CECR1
    CSRP2BP
    DEFB131
    ERBB2IP
    FAM175A
    FAM175B
    FIGF
    FYB
    GYLTL1B
    IKBKAP
    INADL
    KIAA0101
    KIRREL
    MINA
    MRE11A
    NGFRAP1
    PARK2
    PTRF
    RQCD1
    SEPP1
    SHFM1
    TCEB1
    TCEB2
    TCEB3
    TCEB3B
    TCEB3C
    TCEB3CL
    UFD1L
    VPRBP
    WBSCR17
    WHSC1
    WHSC1L1
    '''
}

class SLAPEnrichEvidenceGenerator():
    def __init__(self, inputFile, mappingStep):
        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('SLAPEnrich') \
                .getOrCreate()

        # Initialize mapping variables
        self.mappingStep = mappingStep
    
        # Initialize input files
        self.inputFile = inputFile
        self.dataframe = None

    def writeEvidenceFromSource(self):
        '''
        Processing of the input file to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''
        # Read input file
        self.dataframe = self.spark \
                        .read.csv("slapenrich_opentargets.tsv", sep=r'\t', header=True) \
                        .select("ctype", "gene", "pathway", "SLAPEnrichPval") \
                        .withColumnRenamed("ctype", "Cancer_type_acronym") \
                        .withColumnRenamed("SLAPEnrichPval", "pval")

        # Mapping step
        if self.mappingStep:
            self.dataframe = self.cancer2EFO()

        # Build evidence strings per row
        evidences = self.dataframe.rdd \
            .map(SLAPEnrichEvidenceGenerator.parseEvidenceString) \
            .collect() # list of dictionaries
        
        return evidences
    
    def cancer2EFO(self):
        diseaseMappingsFile = self.spark \
                        .read.csv("resources/cancer2EFO_mappings.tsv", sep=r'\t', header=True) \
                        .select("Cancer_type_acronym", "EFO_id") \

        self.dataframe = self.dataframe.join(
            diseaseMappingsFile,
            on="Cancer_type_acronym",
            how="inner"
        )

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "slapenrich",
                "datatypeId" : "affected_pathway",
                "diseaseFromSource" : row["Cancer_type_acronym"],
                "diseaseFromSourceMappedId" : row["EFO_id"],
                "resourceScore" : row["pval"],
                "pathwayName" : row["pathway"], #split pathway
                "pathwayId" : row["pathway"], #split pathway
                "targetFromSourceId" : row["gene"]
            }
            return evidence
        except Exception as e:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences for the SLAPEnrich data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input source .tsv file.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the evidence compressed JSON file containing the evidence strings.")
    parser.add_argument("-m", "--mappingStep", required=False, type=bool, default=True, help="State whether to run the disease to EFO term mapping step.")

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    outputFile = args.outputFile
    mappingStep = args.mappingStep

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidenceBuilder = SLAPEnrichEvidenceGenerator(inputFile, mappingStep)

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.writeEvidenceFromSource()

    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()