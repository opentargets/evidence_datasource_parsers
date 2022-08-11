#!/usr/bin/env python

import argparse
import logging
import sys

from pyspark.sql.types import StringType
import pyspark.sql.functions as F

from common.evidence import initialize_sparksession, read_path, write_evidence_strings


# The following target labels are excluded as they were grounded to too many target Ids
EXCLUDED_TARGET_TERMS = [
    'TEC',
    'TECS',
    'Tec',
    'tec',
    '\'',
    '(',
    ')',
    '-',
    '-S',
    'S',
    'S-',
    'SS',
    'SSS',
    'Ss',
    'Ss-',
    's',
    's-',
    'ss',
    'U3',
    'U6',
    'u6',
    'SNORA70',
    'U2',
    'U8',
]

# The following are the relevant sections for disease/target associations as described in PMID28587637
SECTIONS_OF_INTEREST = [
    "title",
    "abstract",
    "intro",
    "case",
    "figure",
    "table",
    "discuss",
    "concl",
    "results",
    "appendix",
    "other",
]


def main(cooccurrences, outputFile):

    # Log parameters:
    logging.info(f'Cooccurrence file: {cooccurrences}')
    logging.info(f'Output file: {outputFile}')
    logging.info('Generating evidence:')

    # Load/filter datasets:
    agg_cooccurrence_df = (
        # Reading file:
        read_path(cooccurrences, spark)
        .repartition(200)
        # Filter out pairs found in unwanted sections
        .filter(F.col('section').isin(SECTIONS_OF_INTEREST))
        # Casting integer pmid column to string:
        .withColumn("pmid", F.trim(F.col('pmid').cast(StringType())))
        # Dropping pmcid values that violate schema:
        .withColumn('pmcid', F.when(F.col('pmcid').rlike(r'^PMC\d+$'), F.col('pmcid')))
        # Publication identifier is a pmid if available, otherwise pmcid
        .withColumn('publicationIdentifier', F.when(F.col('pmid').isNull(), F.col('pmcid')).otherwise(F.col('pmid')))
        # Filtering for disease/target cooccurrences:
        .filter(
            (F.col('type') == 'GP-DS')  # Filter gene/protein - disease cooccurrence
            & F.col('isMapped')  # Filtering for mapped cooccurrences
            & F.col('publicationIdentifier').isNotNull() # Making sure at least the pmid or the pmcid is given:
            & (F.length(F.col('text')) < 600) # Exclude sentences with more than 600 characters
            & (F.col('label1').isin(EXCLUDED_TARGET_TERMS) == False)  # Excluding target labels from the exclusion list
        )
        # Renaming columns:
        .withColumnRenamed('keywordId1', 'targetFromSourceId')
        .withColumnRenamed('keywordId2', 'diseaseFromSourceMappedId')
        # Aggregating data by publication, target and disease:
        .groupBy(['publicationIdentifier', 'targetFromSourceId', 'diseaseFromSourceMappedId'])
        .agg(
            F.collect_set(F.col('pmcid')).alias('pmcIds'),
            F.collect_set(F.col('pmid')).alias('literature'),
            F.min('year').alias('publicationYear'),
            F.collect_set(
                F.struct(
                    F.col('text'),
                    F.col('start1').alias('tStart'),
                    F.col('end1').alias('tEnd'),
                    F.col('start2').alias('dStart'),
                    F.col('end2').alias('dEnd'),
                    F.col('section'),
                )
            ).alias('textMiningSentences'),
            F.sum(F.col('evidence_score')).alias('resourceScore'),
        )
        # Nullify pmcIds if empty array:
        .withColumn('pmcIds', F.when(F.size('pmcIds') != 0, F.col('pmcIds')))
        # Only evidence with score above 1 is considered:
        .filter(F.col('resourceScore') > 1)
    )

    # Final formatting and saving data:
    evidence = (
        agg_cooccurrence_df
        # Adding literal columns:
        .withColumn('datasourceId', F.lit('europepmc')).withColumn('datatypeId', F.lit('literature'))
        # Reorder columns:
        .select(
            [
                'datasourceId',
                'datatypeId',
                'targetFromSourceId',
                'diseaseFromSourceMappedId',
                'resourceScore',
                'literature',
                'textMiningSentences',
                'pmcIds',
                'publicationYear',
            ]
        )
    )

    write_evidence_strings(evidence, outputFile)
    logging.info('EPMC disease target evidence saved.')
    logging.info(f'Number of evidence: {agg_cooccurrence_df.count()}')
    # Report on the number of diseases, targets and associations if loglevel == "debug" to avoid cost on computation time:
    logging.debug(f"Number of publications: {agg_cooccurrence_df.select(F.col('publicationIdentifier')).count()}")
    logging.debug(
        f"Number of publications without pubmed ID: {agg_cooccurrence_df.filter(F.col('publicationIdentifier').contains('PMC')).select('publicationIdentifier').distinct().count()}"
    )
    logging.debug(f"Number of targets: {evidence.select(F.col('targetFromSourceId')).distinct().count()}")
    logging.debug(f"Number of diseases: {evidence.select(F.col('diseaseFromSourceMappedId')).distinct().count()}")
    logging.debug(
        f"Number of associations: {evidence.select(F.col('diseaseFromSourceMappedId'), F.col('targetFromSourceId')).dropDuplicates().count()}"
    )


def parse_args():

    parser = argparse.ArgumentParser(
        description='This script generates target/disease evidence strings from ePMC cooccurrence files.'
    )
    parser.add_argument(
        '--cooccurrences', help='Directory of parquet with the ePMC cooccurrences', type=str, required=True
    )
    parser.add_argument(
        '--outputFile', help='Resulting evidence file saved as compressed JSON.', type=str, required=True
    )
    parser.add_argument('--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    args = parser.parse_args()

    # extract parameters:
    cooccurrences = args.cooccurrences
    logFile = args.logFile
    outputFile = args.outputFile

    return (cooccurrences, logFile, outputFile)


if __name__ == '__main__':

    # Parse arguments:
    cooccurrences, logFile, outputFile = parse_args()

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

    global spark
    spark = initialize_sparksession()

    # Calling main function:
    main(cooccurrences, outputFile)
