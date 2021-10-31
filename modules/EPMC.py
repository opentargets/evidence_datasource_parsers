#!/usr/bin/env python

import argparse
import logging
import sys

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.types import StringType
import pyspark.sql.functions as pf

from common.spark import detect_spark_memory_limit


# The following target labels are excluded as they were grounded to too many target Ids
EXCLUDED_TARGET_TERMS = ['TEC', 'TECS', 'Tec', 'tec', '\'', '(', ')', '-', '-S', 'S', 'S-', 'SS', 'SSS',
    'Ss', 'Ss-', 's', 's-', 'ss', 'U3', 'U6', 'u6', 'SNORA70', 'U2', 'U8']


def main(cooccurrenceFile, outputFile, local=False):

    # Initialize spark session
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
        .master('local[*]')
        .getOrCreate()
    )

    logging.info(f'Spark version: {spark.version}')

    # Log parameters:
    logging.info(f'Cooccurrence file: {cooccurrenceFile}')
    logging.info(f'Output file: {outputFile}')
    logging.info('Generating evidence:')

    # Load/filter datasets:
    filtered_cooccurrence_df = (
        # Reading file:
        spark.read.parquet(cooccurrenceFile)

        # Casting integer pmid column to string:
        .withColumn("pmid", pf.col('pmid').cast(StringType()))

        # Publication identifier is a pmid if available, otherwise pmcid
        .withColumn(
            'publicationIdentifier',
            pf.when(pf.col('pmid').isNull(), pf.col('pmcid'))
            .otherwise(pf.col('pmid'))
        )

        # Filtering for disease/target cooccurrences:
        .filter(
            (pf.col('type') == 'GP-DS') &  # Filter gene/protein - disease cooccurrence
            (pf.col('isMapped')) &  # Filtering for mapped cooccurrences
            (pf.col('publicationIdentifier').isNotNull()) &  # Making sure at least the pmid or the pmcid is given:
            (pf.length(pf.col('text')) < 600) &  # Exclude sentences with more than 600 characters
            (pf.col('label1').isin(EXCLUDED_TARGET_TERMS) == False)  # Excluding target labels from the exclusion list
        )

        # Renaming columns:
        .withColumnRenamed('keywordId1', 'targetFromSourceId')
        .withColumnRenamed('keywordId2', 'diseaseFromSourceMappedId')
    )

    # Report on the number of diseases, targets and associations if loglevel == "debug" to avoid cost on computation time:
    logging.debug(f"Number of publications: {filtered_cooccurrence_df.select(pf.col('publicationIdentifier')).distinct().count()}")
    logging.debug(f"Number of targets: {filtered_cooccurrence_df.select(pf.col('targetFromSourceId')).distinct().count()}")
    logging.debug(f"Number of diseases: {filtered_cooccurrence_df.select(pf.col('diseaseFromSourceMappedId')).distinct().count()}")
    logging.debug(f"Number of associations: {filtered_cooccurrence_df.select(pf.col('diseaseFromSourceMappedId'), pf.col('targetFromSourceId')).dropDuplicates().count()}")
    logging.debug(f"Number of publications without pubmed ID: {filtered_cooccurrence_df.filter(pf.col('pmid').isNull()).select('pmcid').distinct().count()}")

    # Aggregating cooccurrence, get score apply filter:
    aggregated_df = (
        filtered_cooccurrence_df

        # Aggregating data by publication, target and disease:
        .groupBy(['publicationIdentifier', 'targetFromSourceId', 'diseaseFromSourceMappedId'])
        .agg(
            pf.collect_set(pf.col('pmcid')).alias('pmcIds'),
            pf.collect_set(pf.col('pmid')).alias('literature'),
            pf.collect_set(
                pf.struct(
                    pf.col('text'),
                    pf.col('start1').alias('tStart'),
                    pf.col('end1').alias('tEnd'),
                    pf.col('start2').alias('dStart'),
                    pf.col('end2').alias('dEnd'),
                    pf.col('section')
                )
            ).alias('textMiningSentences'),
            pf.sum(pf.col('evidence_score')).alias('resourceScore')
        )

        # Only evidence with score above 1 is considered:
        .filter(pf.col('resourceScore') > 1)
    )

    # Report number of evidence:
    logging.info(f'Number of evidence: {aggregated_df.count()}')

    # Final formatting and saving data:
    (
        aggregated_df

        # Adding literal columns:
        .withColumn('datasourceId', pf.lit('europepmc'))
        .withColumn('datatypeId', pf.lit('literature'))

        # Reorder columns:
        .select(['datasourceId', 'datatypeId', 'targetFromSourceId', 'diseaseFromSourceMappedId', 'resourceScore',
                 'literature', 'textMiningSentences', 'pmcIds'])

        # Save output:
        .write.format('json').mode('overwrite').option('compression', 'gzip').save(outputFile)
    )

    logging.info('EPMC disease target evidence saved.')


def parse_args():

    parser = argparse.ArgumentParser(
        description='This script generates target/disease evidence strings from ePMC cooccurrence files.')
    parser.add_argument(
        '--cooccurrenceFile', help='Partioned parquet file with the ePMC cooccurrences', type=str, required=True)
    parser.add_argument(
        '--outputFile', help='Resulting evidence file saved as compressed JSON.', type=str, required=True)
    parser.add_argument(
        '--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    parser.add_argument(
        '--local', help='Destination of the logs generated by this script.', action='store_true', required=False, default=False)
    args = parser.parse_args()

    # extract parameters:
    cooccurrenceFile = args.cooccurrenceFile
    logFile = args.logFile
    outputFile = args.outputFile
    local = args.local

    return (cooccurrenceFile, logFile, outputFile, local)


if __name__ == '__main__':

    # Parse arguments:
    cooccurrenceFile, logFile, outputFile, local = parse_args()

    # Initialize logger based on the provided logfile.
    # If no logfile is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if logFile:
        logging.config.fileConfig(filename=logFile)
    else:
        logging.StreamHandler(sys.stderr)

    # Calling main function:
    main(cooccurrenceFile, outputFile, local)
