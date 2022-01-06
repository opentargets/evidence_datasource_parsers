#!/usr/bin/env python3
"""Parser for Target Enabling Packages (TEP) downloaded from the Structural Genomics Consortium website."""

import argparse
import logging
import sys

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.functions import explode, col, concat, lit, lower, regexp_replace, split
from pyspark import SparkFiles

from common.evidence import detect_spark_memory_limit, write_evidence_strings


'''

# The TEP dataset is made available by the SGC as a tab-separated file:
TEPURL = 'https://www.thesgc.org/sites/default/files/fileuploads/available-teps.tsv'


def main(outputFile: str) -> None:

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
        .config("spark.sql.broadcastTimeout", "36000")
        .master('local[*]')
        .getOrCreate()
    )

    spark.sparkContext.addFile("https://www.thesgc.org/sites/default/files/fileuploads/available-teps.tsv")

    # Fetching and processing the TEP table and saved as a JSON file:
    TEP_df = (
        spark.read.csv(SparkFiles.get(TEPURL.split('/')[-1]), sep='\t', header=True)

        # Generating TEP url from Gene column: SLC12A4/SLC12A6 -> https://www.thesgc.org/tep/SLC12A4SLC12A6
        .withColumn('url', concat(lit('https://www.thesgc.org/tep/'), regexp_replace(lower(col('Gene')), '/', '')))

        # Exploding TEPs, where multiple genes are given:
        .withColumn('targetFromSourceId', explode(split(col('Gene'), '/')))

        # Renaming columns:
        .withColumnRenamed('Therapeutic Area', 'therapeuticArea')
        .withColumnRenamed('Description', 'description')

        # Dropping columns:
        .drop(*['Gene', 'version', 'Date'])
        .persist()
    )

    logging.info('TEP dataset has been downloaded and formatted.')
    logging.info(f'Number of TEPs: {TEP_df.count()}')
    logging.info(f'Number of unique genes: {TEP_df.select("targetFromSourceId").distinct().count()}')

    # Saving data:
    write_evidence_strings(TEP_df, outputFile)
    logging.info(f'TEP dataset is written to {outputFile}.')


if __name__ == '__main__':

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(description='This script fetches TEP data from Structural Genomics Consortium.')
    parser.add_argument('--output_file', '-o', type=str, help='Output file. gzipped JSON', required=True)
    parser.add_argument('--log_file', type=str, help='File into which the logs are saved', required=False)
    args = parser.parse_args()

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    main(args.output_file)
