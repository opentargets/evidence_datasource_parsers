#!/usr/bin/env python3
"""Evidence parser for Gene2Phenotype's disease panels."""

import argparse
import logging
import os
import sys
import tempfile

from ontoma import OnToma
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import array, col, concat, explode, first, lit, split, udf, when
from pyspark.sql.types import ArrayType, IntegerType, StringType, StructField, StructType, TimestampType

from common.ontology import add_efo_mapping

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

class ontoma_efo_lookup():
    """
    Map orphanet diseases to the EFO ontology
    """
    def __init__(self):
        self.otmap = OnToma(cache_dir='cache')

    def get_mapping(self, terms=[]):
        disease_label, disease_id = terms
        label_mapping = [
            result.id_ot_schema
            for result in self.otmap.find_term(disease_label)]

        if len(label_mapping) != 0:
            return (disease_label, disease_id, label_mapping)
        try:
            id_mapping = [
                result.id_ot_schema
                for result in self.otmap.find_term(disease_id)]
            return (disease_label, disease_id, id_mapping)
        except AttributeError:
            return None

def main(dd_file, eye_file, skin_file, cancer_file, output_file, local):

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

    # Initialize disease mapping object:
    ol_obj = ontoma_efo_lookup()

    # Specify schema -> this schema is applied for all gene2phenotype files:
    gene2phenotype_schema = (
        StructType()
        .add('gene symbol', StringType())
        .add('gene mim', IntegerType())
        .add('disease name', StringType())
        .add('disease mim', StringType())
        .add('DDD category', StringType())
        .add('allelic requirement', StringType())
        .add('mutation consequence', StringType())
        .add('phenotypes', StringType())
        .add('organ specificity list', StringType())
        .add('pmids', StringType())
        .add('panel', StringType())
        .add('prev symbols', StringType())
        .add('hgnc id', IntegerType())
        .add('gene disease pair entry date', TimestampType())
    )

    # Load all files for one go:
    gene2phenotype_data = (
        spark.read.csv(
            [dd_file, eye_file, skin_file, cancer_file], schema=gene2phenotype_schema, enforceSchema=True, header=True
        )
        # Split pubmed IDs to list:
        .withColumn('literature', split(col('pmids'), ';'))

        # Split phenotypes:
        .withColumn('phenotypes', split(col('phenotypes'), ';'))

        # Split organ specificity:
        .withColumn('organ_specificities', split(col('organ specificity list'), ';'))
    )

    # Processing data:
    evidence_df = (
        gene2phenotype_data

        # Renaming columns:
        .withColumnRenamed('gene symbol', 'targetFromSourceId')
        .withColumnRenamed('disease name', 'diseaseFromSource')
        .withColumnRenamed('panel', 'studyId')
        .withColumnRenamed('DDD category', 'confidence')

        .withColumn(
            'allelicRequirements',
            when(
                col('allelic requirement').isNotNull(),
                array(col('allelic requirement'))
            ))

        # Process diseaseFromSource
        .withColumn(
            'diseaseFromSourceId',
            when(~col('disease mim').contains('No disease mim'), col('disease mim'))
        )
        .withColumn('diseaseFromSourceId', concat(lit('OMIM:'), col('diseaseFromSourceId')))

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
    )

    '''

    # Get all the diseases + map disease to EFO:
    diseases = (
        evidence_df
        .select('diseaseFromSource', 'diseaseFromSourceId')
        .distinct()
        .collect()
    )

    mapped_diseases = [ol_obj.get_mapping(x) for x in diseases]

    schema = StructType([
        StructField('diseaseFromSource', StringType(), True),
        StructField('diseaseFromSourceId', StringType(), True),
        StructField('diseaseFromSourceMappedId', ArrayType(StringType()), True)
    ])
    mapped_diseases_df = (
        # Dataframe from list of label/id/mapped_id tuples
        spark.createDataFrame(mapped_diseases, schema=schema)
        # Coalesce df row-wise
        .groupBy(
            'diseaseFromSource', 'diseaseFromSourceId')
        .agg(first("diseaseFromSourceMappedId", ignorenulls=True).alias("diseaseFromSourceMappedId"))
        # Explode cases where one diseaseFromSource maps to many EFO terms
        .withColumn('diseaseFromSourceMappedId', explode('diseaseFromSourceMappedId'))
    )

    # Merge evidence with the mapped disease:
    evidence_df = (
        evidence_df
        .join(mapped_diseases_df, on=['diseaseFromSource', 'diseaseFromSourceId'], how='left')
    )
    '''

    evidence_df = add_efo_mapping(evidence_strings=evidence_df, spark_instance=spark)

    # Saving data:
    write_evidence_strings(evidence_df, output_file)
    logging.info(f'{evidence_df.count()} evidence strings have been saved to {output_file}')


def write_evidence_strings(evidence: DataFrame, output_file: str) -> None:
    '''Exports the table to a compressed JSON file containing the evidence strings'''
    with tempfile.TemporaryDirectory() as tmp_dir_name:
        (
            evidence.coalesce(1).write.format('json').mode('overwrite')
            .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
        )
        json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
        assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
        os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-d', '--dd_panel',
                        help='DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-e', '--eye_panel',
                        help='Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-s', '--skin_panel',
                        help='Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('-c', '--cancer_panel',
                        help='Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        type=str)
    parser.add_argument('--local', help='Where the ', action='store_true', required=False, default=False)
    parser.add_argument('-o', '--output_file', help='Absolute path of the gzipped, JSON evidence file.', type=str)
    parser.add_argument('-l', '--log_file', help='Filename to store the parser logs.', type=str)

    args = parser.parse_args()

    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    output_file = args.output_file
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
    main(dd_file, eye_file, skin_file, cancer_file, output_file, local)
