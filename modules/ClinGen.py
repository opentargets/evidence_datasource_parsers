#!/usr/bin/env python3
"""Evidence parser for ClinGen's Gene Validity Curations."""

import argparse
import logging
import os
import tempfile

from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import array, col, first, lit, monotonically_increasing_id, struct, trim
from pyspark.sql.types import StringType, StructType, TimestampType

from common.ontology import add_efo_mapping
from common.spark import write_evidence_strings


def main(input_file: str, output_file: str, cache_dir: str, local: bool = False) -> None:

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
            .getOrCreate()
        )

    # Read and process Clingen's table into evidence strings

    clingen_df = read_input_file(input_file, spark_instance=spark)
    logging.info('Gene Validity Curations table has been imported. Processing evidence strings.')

    evidence_df = process_clingen(clingen_df)

    evidence_df = add_efo_mapping(evidence_strings=evidence_df, spark_instance=spark, ontoma_cache_dir=cache_dir)
    logging.info('Disease mappings have been added.')

    write_evidence_strings(evidence_df, output_file)
    logging.info(f'{evidence_df.count()} evidence strings have been saved to {output_file}')


def read_input_file(input_file: str, spark_instance) -> DataFrame:
    """
    Reads Gene Validity Curations CSV file into a Spark DataFrame forcing the schema
    The first 6 rows of this file include metadata that needs to be dropped
    """

    clingen_schema = (
        StructType()
        .add('GENE SYMBOL', StringType())
        .add('GENE ID (HGNC)', StringType())
        .add('DISEASE LABEL', StringType())
        .add('DISEASE ID (MONDO)', StringType())
        .add('MOI', StringType())
        .add('SOP', StringType())
        .add('CLASSIFICATION', StringType())
        .add('ONLINE REPORT', StringType())
        .add('CLASSIFICATION DATE', TimestampType())
        .add('GCEP', StringType())
    )

    clingen_df = (
        spark_instance.read.csv(input_file, schema=clingen_schema)
        # The first 6 rows of the GVC file include metadata that needs to be dropped
        .withColumn('idx', monotonically_increasing_id())
        .filter(col('idx') > 5)
        .drop('idx')
    )

    return clingen_df


def process_clingen(clingen_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """

    return (
        clingen_df
        .withColumn('datasourceId', lit('clingen'))
        .withColumn('datatypeId', lit('genetic_literature'))
        .withColumn('urls', struct(col('ONLINE REPORT').alias('url')))

        .select(
            'datasourceId', 'datatypeId',
            trim(col('GENE SYMBOL')).alias('targetFromSourceId'),
            col('DISEASE LABEL').alias('diseaseFromSource'),
            col('DISEASE ID (MONDO)').alias('diseaseFromSourceId'),
            array(col('MOI')).alias('allelicRequirements'),
            array(col('urls')).alias('urls'),
            col('CLASSIFICATION').alias('confidence'),
            col('GCEP').alias('studyId')
        )
    )


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse ClinGen gene-disease associations from Gene Validity Curations')
    parser.add_argument('--input_file',
                        help='Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity',
                        type=str, required=True)
    parser.add_argument('--output_file',
                        help='Absolute path of the gzipped JSON evidence file.',
                        type=str, required=True)
    parser.add_argument('--log_file', type=str,
                        help='Optional filename to redirect the logs into.')
    parser.add_argument('--cache_dir', required=False, help='Directory to store the OnToma cache files in.')
    parser.add_argument(
        '--local', action='store_true', required=False, default=False,
        help='Flag to indicate if the script is executed locally or on the cluster'
    )
    args = parser.parse_args()

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)

    # Report input data:
    logging.info(f'Clingen input file path: {args.input_file}')
    logging.info(f'Evidence output file path: {args.output_file}')
    logging.info(f'Cache directory: {args.cache_dir}')

    main(
        input_file=args.input_file, output_file=args.output_file,
        cache_dir=args.cache_dir, local=args.local
    )
