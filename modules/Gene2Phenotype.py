#!/usr/bin/env python3
"""Evidence parser for Gene2Phenotype's disease panels."""

import argparse
import logging
import sys

from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import array, col, concat, lit, split, udf, when
from pyspark.sql.types import IntegerType, StringType, StructType, TimestampType

from common.ontology import add_efo_mapping
from common.evidence import detect_spark_memory_limit, write_evidence_strings


G2P_mutationCsq2functionalCsq = {
    'uncertain': 'SO_0002220',  # function_uncertain_variant
    'absent gene product': 'SO_0002317',  # absent_gene_product
    'altered gene product structure': 'SO_0002318',  # altered_gene_product_structure
    '5_prime or 3_prime UTR mutation': 'SO_0001622',  # UTR_variant
    'increased gene product level': 'SO_0002315',  # increased_gene_product_level
    'cis-regulatory or promotor mutation': 'SO_0001566',  # regulatory_region_variant
}


def main(
    dd_file: str, eye_file: str, skin_file: str, cancer_file: str, cardiac_file: str, output_file: str, cache_dir: str, local: bool = False
) -> None:

    # Initialize spark session
    global spark
    spark_mem_limit = detect_spark_memory_limit()
    spark_conf = (
        SparkConf()
        .set('spark.driver.memory', f'{spark_mem_limit}g')
        .set('spark.executor.memory', f'{spark_mem_limit}g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
    )
    spark = (
        SparkSession.builder
        .config(conf=spark_conf)
        .config("spark.sql.broadcastTimeout", "36000")
        .master('local[*]')
        .getOrCreate()
    )

    # Read and process G2P's tables into evidence strings
    gene2phenotype_df = read_input_file(
        dd_file, eye_file, skin_file, cancer_file, cardiac_file)
    logging.info('Gene2Phenotype panels have been imported. Processing evidence strings.')

    evidence_df = process_gene2phenotype(gene2phenotype_df)

    evidence_df = add_efo_mapping(evidence_strings=evidence_df, ontoma_cache_dir=cache_dir, spark_instance=spark)
    logging.info('Disease mappings have been added.')

    # Saving data:
    write_evidence_strings(evidence_df, output_file)
    logging.info(f'{evidence_df.count()} evidence strings have been saved to {output_file}')

def read_input_file(
    dd_file: str, eye_file: str, skin_file: str, cancer_file: str, cardiac_file: str
) -> DataFrame:
    '''
    Reads G2P's panel CSV files into a Spark DataFrame forcing the schema
    '''

    gene2phenotype_schema = (
        StructType()
        .add('gene symbol', StringType())
        .add('gene mim', IntegerType())
        .add('disease name', StringType())
        .add('disease mim', StringType())
        .add('confidence category', StringType())
        .add('allelic requirement', StringType())
        .add('mutation consequence', StringType())
        .add('phenotypes', StringType())
        .add('organ specificity list', StringType())
        .add('pmids', StringType())
        .add('panel', StringType())
        .add('prev symbols', StringType())
        .add('hgnc id', IntegerType())
        .add('gene disease pair entry date', TimestampType())
        .add('cross cutting modifier', StringType())
        .add('mutation consequence flag', StringType())
    )

    return (
        spark.read.csv(
            [dd_file, eye_file, skin_file, cancer_file, cardiac_file], schema=gene2phenotype_schema, enforceSchema=True, header=True
        )
    )

def process_gene2phenotype(gene2phenotype_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """

    evidence_df = (
        gene2phenotype_df

        # Split pubmed IDs to list:
        .withColumn('literature', split(col('pmids'), ';'))

        # Split phenotypes:
        .withColumn('phenotypes', split(col('phenotypes'), ';'))

        # Split organ specificity:
        .withColumn('organ_specificities', split(col('organ specificity list'), ';'))

        # Reshaping columns:
        .withColumnRenamed('gene symbol', 'targetFromSourceId')
        .withColumnRenamed('disease name', 'diseaseFromSource')
        .withColumnRenamed('panel', 'studyId')
        .withColumnRenamed('confidence category', 'confidence')
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

    return evidence_df

def translate(mapping):
    """
    Mapping consequences - to SO codes
    """
    def translate_(col):
        return mapping.get(col)
    return udf(translate_, StringType())


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description='Parse Gene2Phenotype gene-disease files downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/')
    parser.add_argument('-d', '--dd_panel',
                        help='DD panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        required=True, type=str)
    parser.add_argument('-e', '--eye_panel',
                        help='Eye panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        required=True, type=str)
    parser.add_argument('-s', '--skin_panel',
                        help='Skin panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        required=True, type=str)
    parser.add_argument('-c', '--cancer_panel',
                        help='Cancer panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        required=True, type=str)
    parser.add_argument('-cr', '--cardiac_panel',
                        help='Cardiac panel file downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads',
                        required=True, type=str)
    parser.add_argument('-o', '--output_file', help='Absolute path of the gzipped, JSON evidence file.', required=True, type=str)
    parser.add_argument('-l', '--log_file', help='Filename to store the parser logs.', type=str)
    parser.add_argument('--cache_dir', required=False, help='Directory to store the OnToma cache files in.')
    parser.add_argument(
        '--local', action='store_true', required=False, default=False,
        help='Flag to indicate if the script is executed locally or on the cluster'
    )

    args = parser.parse_args()

    # Get parameters
    dd_file = args.dd_panel
    eye_file = args.eye_panel
    skin_file = args.skin_panel
    cancer_file = args.cancer_panel
    cardiac_file = args.cardiac_panel
    output_file = args.output_file
    log_file = args.log_file
    cache_dir = args.cache_dir
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
    logging.info(f'Cardiac panel file: {cardiac_file}')

    # Calling main:
    main(dd_file, eye_file, skin_file, cancer_file, cardiac_file, output_file, cache_dir, local)
