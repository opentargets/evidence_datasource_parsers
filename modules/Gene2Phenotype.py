import logging
import csv
import gzip
import argparse
import json
import sys

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.functions import split, col, udf, lit
from pyspark.sql.types import StringType, IntegerType, TimestampType, StructType

import ontoma


G2P_mutationCsq2functionalCsq = {
    'loss of function': 'SO_0002054',  # loss_of_function_variant
    'all missense/in frame': 'SO_0001650',  # inframe_variant
    'uncertain': 'SO_0002220',  # function_uncertain_variant
    'activating': 'SO_0002053',  # gain_of_function_variant
    'dominant negative': 'SO_0002052',  # dominant_negative_variant
    '': None,
    'gain of function': 'SO_0002053',  # gain_of_function_variant
    'cis-regulatory or promotor mutation': 'SO_0001566',  # regulatory_region_variant
    '5_prime or 3_prime UTR mutation': 'SO_0001622',  # UTR_variant
    'increased gene dosage': 'SO_0001911',  # copy_number_increase
    'part of contiguous gene duplication': 'SO_1000173'  # tandem_duplication
}


def translate(mapping):
    '''
    Mapping consequences - to SO codes
    '''
    def translate_(col):
        return mapping.get(col)
    return udf(translate_, StringType())

class disease_map(object):
    
    def __init__(self):
        self.ontoma = ontoma.interface.OnToma()

    def map_disease(self, disease_name, omim_id):
        logging.info(f"Mapping '{disease_name}'")

        # Search disease name using OnToma and accept perfect matches
        ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
        
        # If there's some mapping available:
        if ontoma_mapping:
            
            # Extracting term if no action is required:
            if ontoma_mapping['action'] is None:
                return ontoma_mapping
                
            # When there is an exact match, but action is required:
            elif ontoma_mapping['quality'] == "match":
                
                # Match in HP or ORDO, check if there is a match in MONDO too. If so, give preference to MONDO hit
                mondo_mapping = self.search_mondo(disease_name)
                
                if mondo_mapping:
                    # Mondo mapping good - return
                    if mondo_mapping['exact']:
                        return mondo_mapping
                    # Mondo mapping bad - return ontoma
                    else:
                        return ontoma_mapping 
                else:
                    # Mondo mapping bad - return ontoma
                    return ontoma_mapping

            else:
                # OnToma fuzzy match. First check if the mapping term has a xref to the OMIM id. 
                # If not, check in MONDO and if there is not match ignore evidence and report disease
                if self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                    for efo_xref in self.ontoma.get_efo_from_xref(f"OMIM:{omim_id}"):
                        # Extract EFO id from OnToma results
                        efo_id = ontoma_mapping['term'].split('/')[-1].replace('_', ':')

                        if efo_id == efo_xref['id']:
                            return ontoma_mapping

                # xref search didn't work, try MONDO as the last resort
                mondo_mapping = self.search_mondo(disease_name)
                if mondo_mapping:
                    if mondo_mapping['exact']:
                        return mondo_mapping
                    else:
                        return None
                else:
                    # Record the unmapped disease
                    return None

        else:
            # No match in EFO, HP or ORDO
            mondo_mapping = self.search_mondo(disease_name)
            if mondo_mapping:
                if mondo_mapping['exact']:
                    return mondo_mapping
                else:
                    return None
            else:
                return None

            
    def search_mondo(self, disease_name):

        disease_name = disease_name.lower()

        # mondo_lookup works like a dictionary lookup so if disease is not in there it raises and error instead of returning `None`
        try:
            mondo_term = self.ontoma.mondo_lookup(disease_name)
            return {
                'id': mondo_term, 
                'name': self.ontoma.get_mondo_label(mondo_term), 
                'exact': True
            }
        except KeyError as e:
            exact_ols_mondo = self.ontoma._ols.besthit(disease_name, ontology=['mondo'], field_list=['iri', 'label'], exact=True)
            
            if exact_ols_mondo:
                return {'term': exact_ols_mondo['iri'], 'name': exact_ols_mondo['label'], 'exact':True}
            
            else:
                ols_mondo = self.ontoma._ols.besthit(disease_name,
                                                     ontology=['mondo'],
                                                     field_list=['iri', 'label'],
                                                     bytype='class')
                if ols_mondo:
                    return {'term': ols_mondo['iri'], 'name': ols_mondo['label'], 'exact': False}
                else:
                    return None

def main(dd_file, eye_file, skin_file, cancer_file, outfile, local):

    # Initialize disease mapping object:
    dm_obj = disease_map()

    # UDF to look up EFO mappings:
    @udf(StringType())
    def map_disease(label, disease_id):
        lookup = dm_obj.map_disease(label, disease_id)
        if lookup:
            try:
                return lookup['term'].split('/')[-1]
            except:
                print(lookup)
        else:
            return None

    # Initialize spark session
    if local:
        sparkConf = (
            SparkConf()
            .set('spark.driver.memory', '15g')
            .set('spark.executor.memory', '15g')
            .set('spark.driver.maxResultSize', '0')
            .set('spark.debug.maxToStringFields', '2000')
            .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        )
        spark = (
            SparkSession.builder
            .config(conf=sparkConf)
            .config("spark.sql.broadcastTimeout", "36000")
            .master('local[*]')
            .getOrCreate()
        )
    else:
        sparkConf = (
            SparkConf()
            .set('spark.driver.maxResultSize', '0')
            .set('spark.debug.maxToStringFields', '2000')
            .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        )
        spark = (
            SparkSession.builder
            .config(conf=sparkConf)
            .config("spark.sql.broadcastTimeout", "36000")
            .getOrCreate()
        )

    # Specify schema -> this schema is applied for all INTOGen files:
    intogen_schema = (
        StructType()
        .add('gene symbol', StringType())
        .add('gene mim', IntegerType())
        .add('disease name', StringType())
        .add('disease mim', StringType())
        .add('DDD category', StringType())
        .add('allelic requirement list', StringType())
        .add('mutation consequence', StringType())
        .add('phenotype list', StringType())
        .add('organ specificity list', StringType())
        .add('pmid list', StringType())
        .add('panel', StringType())
        .add('prev symbol list', StringType())
        .add('hgnc id', IntegerType())
        .add('gene disease entry date', TimestampType())
    )

    # Load all files for one go:
    intogen_data = (
        spark.read.csv(
            [dd_file, eye_file, skin_file, cancer_file], schema=intogen_schema, enforceSchema=True, header=True
        )

        # Split pubmed IDs to list:
        .withColumn('literature', split(col('pmid list'), ';'))

        # Split phenotypes:
        .withColumn('phenotypes', split(col('phenotype list'), ';'))

        # Split organ specificity:
        .withColumn('organ_specificities', split(col('organ specificity list'), ';'))

        # Split allelic requirements:
        .withColumn('allelicRequirements', split(col('allelic requirement list'), ';'))
    )

    # Processing data:
    evidence_df = (
        intogen_data

        # Renaming columns:
        .withColumnRenamed('gene symbol', 'targetFromSourceId')
        .withColumnRenamed('disease mim', 'diseaseFromSourceId')
        .withColumnRenamed('disease name', 'diseaseFromSource')
        .withColumnRenamed('panel', 'studyId')
        .withColumnRenamed('DDD category', 'confidence')

        # Map functional consequences:
        .withColumn("variantFunctionalConsequenceId", translate(G2P_mutationCsq2functionalCsq)("mutation consequence"))

        # Adding literature columns:
        .withColumn('datasourceId', lit('gene2phenotype'))
        .withColumn('datatypeId', lit('genetic_literature'))

        # Selecting relevant columns:
        .select(
            'datasourceId', 'datatypeId', 'targetFromSourceId', 'diseaseFromSource',
            'diseaseFromSourceId', 'confidence', 'studyId', 'literature',
            'allelicRequirements', 'variantFunctionalConsequenceId'
        )
        .coalesce(1)
    )

    # Get all the diseases + map disease to EFO:
    diseases = (
        evidence_df
        .select('diseaseFromSource', 'diseaseFromSourceId')
        .distinct()
        .withColumn('diseaseFromSourceMappedId', map_disease(col('diseaseFromSource'), col('diseaseFromSourceId')))
        .persist()
    )

    # Merge evidence with the mapped disease:
    evidence_df = (
        evidence_df
        .join(diseases, how='left', on=['diseaseFromSource', 'diseaseFromSourceId'])
    )

    # Saving data:
    (
        evidence_df
        .write.format('json').mode('overwrite').option('compression', 'gzip').save(outfile)
    )


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-d', '--dd_panel',
                        help='DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-e', '--eye_panel',
                        help='Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-k', '--skin_panel',
                        help='Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-c', '--cancer_panel',
                        help='Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('--local', help='Where the ', action='store_true', required=False, default=False)
    parser.add_argument('-o', '--output_file', help='Name of gzipped evidence file', type=str)
    parser.add_argument('-l', '--log_file', help='Name of gzipped evidence file', type=str)

    args = parser.parse_args()

    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    outfile = args.output_file
    log_file = args.log_file
    local = args.local

    # Configure logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Report input data:
    logging.info(f'DD panel file: {dd_file}')
    logging.info(f'Eye panel file: {eye_file}')
    logging.info(f'Skin panel file: {skin_file}')
    logging.info(f'Cancer panel file: {cancer_file}')

    # Calling main:
    main(dd_file, eye_file, skin_file, cancer_file, outfile, local)
