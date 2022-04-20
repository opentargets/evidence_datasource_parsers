#!/usr/bin/env python

import argparse
import logging
import sys

from pyspark.sql.types import StringType
import pyspark.sql.functions as pf

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


def main(cooccurrenceFile, outputFile):

    # Log parameters:
    logging.info(f'Cooccurrence file: {cooccurrenceFile}')
    logging.info(f'Output file: {outputFile}')
    logging.info('Generating evidence:')

    # Load/filter datasets:
    agg_cooccurrence_df = (
        # Reading file:
        read_path(cooccurrenceFile, spark)
        .repartition(200)
        # Filter out pairs found in unwanted sections
        .filter(pf.col('section').isin(SECTIONS_OF_INTEREST))
        # Casting integer pmid column to string:
        .withColumn("pmid", pf.trim(pf.col('pmid').cast(StringType())))
        # Publication identifier is a pmid if available, otherwise pmcid
        .withColumn(
            'publicationIdentifier', pf.when(pf.col('pmid').isNull(), pf.col('pmcid')).otherwise(pf.col('pmid'))
        )
        # Filtering for disease/target cooccurrences:
        .filter(
            (pf.col('type') == 'GP-DS')
            & (pf.col('isMapped'))  # Filter gene/protein - disease cooccurrence
            & (pf.col('publicationIdentifier').isNotNull())  # Filtering for mapped cooccurrences
            & (pf.length(pf.col('text')) < 600)  # Making sure at least the pmid or the pmcid is given:
            & (  # Exclude sentences with more than 600 characters
                pf.col('label1').isin(EXCLUDED_TARGET_TERMS) == False
            )  # Excluding target labels from the exclusion list
        )
        # Renaming columns:
        .withColumnRenamed('keywordId1', 'targetFromSourceId')
        .withColumnRenamed('keywordId2', 'diseaseFromSourceMappedId')
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
                    pf.col('section'),
                )
            ).alias('textMiningSentences'),
            pf.sum(pf.col('evidence_score')).alias('resourceScore'),
        )
        # Nullify pmcIds if empty array
        .withColumn('pmcIds', pf.when(pf.size('pmcIds') != 0, pf.col('pmcIds')))
        # Only evidence with score above 1 is considered:
        .filter(pf.col('resourceScore') > 1)
    )

    # Final formatting and saving data:
    evidence = (
        agg_cooccurrence_df
        # Adding literal columns:
        .withColumn('datasourceId', pf.lit('europepmc')).withColumn('datatypeId', pf.lit('literature'))
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
            ]
        )
    )

    write_evidence_strings(evidence, outputFile)
    logging.info('EPMC disease target evidence saved.')
    logging.info(f'Number of evidence: {agg_cooccurrence_df.count()}')
    # Report on the number of diseases, targets and associations if loglevel == "debug" to avoid cost on computation time:
    logging.debug(f"Number of publications: {agg_cooccurrence_df.select(pf.col('publicationIdentifier')).count()}")
    logging.debug(f"Number of targets: {evidence.select(pf.col('targetFromSourceId')).distinct().count()}")
    logging.debug(f"Number of diseases: {evidence.select(pf.col('diseaseFromSourceMappedId')).distinct().count()}")
    logging.debug(
        f"Number of associations: {evidence.select(pf.col('diseaseFromSourceMappedId'), pf.col('targetFromSourceId')).dropDuplicates().count()}"
    )
    logging.debug(
        f"Number of publications without pubmed ID: {agg_cooccurrence_df.filter(pf.col('pmid').isNull()).select('pmcid').distinct().count()}"
    )


def parse_args():

    parser = argparse.ArgumentParser(
        description='This script generates target/disease evidence strings from ePMC cooccurrence files.'
    )
    parser.add_argument(
        '--cooccurrenceFile', help='Partioned parquet file with the ePMC cooccurrences', type=str, required=True
    )
    parser.add_argument(
        '--outputFile', help='Resulting evidence file saved as compressed JSON.', type=str, required=True
    )
    parser.add_argument('--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    args = parser.parse_args()

    # extract parameters:
    cooccurrenceFile = args.cooccurrenceFile
    logFile = args.logFile
    outputFile = args.outputFile

    return (cooccurrenceFile, logFile, outputFile)


if __name__ == '__main__':

    # Parse arguments:
    cooccurrenceFile, logFile, outputFile = parse_args()

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
    main(cooccurrenceFile, outputFile)
