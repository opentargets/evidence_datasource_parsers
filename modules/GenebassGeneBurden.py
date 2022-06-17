#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the Genebass webservice."""

import argparse
import logging
import sys

from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame
from pyspark.sql.functions import coalesce, col, lit, log10, pow, round, when
from pyspark.sql.types import IntegerType

from common.evidence import initialize_sparksession, read_trait_mappings, write_evidence_strings

METHOD_DESC = {
    'pLoF': 'Burden test carried out with rare pLOF variants.',
    'missense|LC': 'Burden test carried out with rare missense variants including low-confidence pLOF and in-frame indels.',
    'synonymous': 'Burden test carried out with rare synonymous variants.',
    'pLoF|missense|LC': 'Burden test carried out with pLOF or missense variants.',
}


def main(genebass_data: str, spark_instance: SparkSession) -> DataFrame:
    """
    This module extracts and processes target/disease evidence from the raw Genebass Portal.
    """
    logging.info(f'File with the Genebass gene burden results: {genebass_data}')

    # Load data
    genebass_df = (
        spark_instance.read.parquet(genebass_data)
        .filter(col('Pvalue_Burden') <= 6.7e-7)
        .select(
            'gene_id',
            'annotation',
            'n_cases',
            'n_controls',
            'trait_type',
            'phenocode',
            'description',
            'Pvalue_Burden',
            'BETA_Burden',
            'SE_Burden',
        )
        # Bring trait to EFO mappings
        .join(
            read_trait_mappings(spark_instance).withColumnRenamed('PROPERTY_VALUE', 'description'),
            on='description',
            how='left',
        )
        .distinct()
        .repartition(20)
        .persist()
    )

    # Write output
    evd_df = parse_genebass_evidence(genebass_df)

    if evd_df.filter(col('resourceScore') == 0).count() != 0:
        logging.error('There are evidence with a P value of 0.')
        raise AssertionError(
            f"There are {evd_df.filter(col('resourceScore') == 0).count()} evidence with a P value of 0."
        )
    if not 8_000 < evd_df.count() < 9_000:
        logging.error(f'Genebass number of evidence are different from expected: {evd_df.count()}')
        raise AssertionError('Genebass number of evidence are different from expected.')
    logging.info(f'{evd_df.count()} evidence strings have been processed.')

    return evd_df


def parse_genebass_evidence(genebass_df: DataFrame) -> DataFrame:
    """
    Parse Genebass's disease/target evidence.
    Args:
        genebass_df: DataFrame with Genebass's portal data
    Returns:
        evd_df: DataFrame with Genebass's data following the t/d evidence schema.
    """
    to_keep = [
        'datasourceId',
        'datatypeId',
        'targetFromSourceId',
        'diseaseFromSource',
        'diseaseFromSourceId',
        'diseaseFromSourceMappedId',
        'pValueMantissa',
        'pValueExponent',
        'beta',
        'betaConfidenceIntervalLower',
        'betaConfidenceIntervalUpper',
        'oddsRatio',
        'oddsRatioConfidenceIntervalLower',
        'oddsRatioConfidenceIntervalUpper',
        'resourceScore',
        'ancestry',
        'ancestryId',
        'projectId',
        'cohortId',
        'studySampleSize',
        'studyCases',
        'statisticalMethod',
        'statisticalMethodOverview',
    ]

    # WARNING: There are some associations with a p-value of 0.0 in Genebass.
    # This is a bug we still have to ellucidate and it might be due to a float overflow.
    # These evidence need to be manually corrected in order not to lose them and for them to pass validation
    # As an interim solution, their p value will equal to the minimum in the evidence set.
    logging.warning(
        f"There are {genebass_df.filter(col('Pvalue_Burden') == 0.0).count()} evidence with a p-value of 0.0."
    )
    minimum_pvalue = (
        genebass_df.filter(col('Pvalue_Burden') > 0.0).agg({'Pvalue_Burden': 'min'}).collect()[0]['min(Pvalue_Burden)']
    )
    genebass_df = genebass_df.withColumn(
        'Pvalue_Burden', when(col('Pvalue_Burden') == 0.0, lit(minimum_pvalue)).otherwise(col('Pvalue_Burden'))
    )

    return (
        genebass_df.withColumn('datasourceId', lit('gene_burden'))
        .withColumn('datatypeId', lit('genetic_association'))
        .withColumn('projectId', lit('Genebass'))
        .withColumn('cohortId', lit('UK Biobank 450k'))
        .withColumn('ancestry', lit('EUR'))
        .withColumn('ancestryId', lit('HANCESTRO_0009'))
        .withColumnRenamed('gene_id', 'targetFromSourceId')
        .withColumnRenamed('description', 'diseaseFromSource')
        .withColumnRenamed('phenocode', 'diseaseFromSourceId')
        .withColumnRenamed('SEMANTIC_TAG', 'diseaseFromSourceMappedId')
        .withColumnRenamed('Pvalue_Burden', 'resourceScore')
        .withColumn('pValueExponent', log10(col('resourceScore')).cast(IntegerType()) - lit(1))
        .withColumn('pValueMantissa', round(col('resourceScore') / pow(lit(10), col('pValueExponent')), 3))
        # Stats are split taking into consideration the type of the trait
        # Those that are not continuous or categorical were reviewed and all of them are considered as categorical
        .withColumn(
            'beta',
            when(col('trait_type') == 'continuous', col('BETA_Burden')),
        )
        .withColumn(
            'betaConfidenceIntervalLower',
            when(col('trait_type') == 'continuous', col('BETA_Burden') - col('SE_Burden')),
        )
        .withColumn(
            'betaConfidenceIntervalUpper',
            when(col('trait_type') == 'continuous', col('BETA_Burden') + col('SE_Burden')),
        )
        .withColumn(
            'oddsRatio',
            when(col('trait_type').isin(['categorical', 'icd_first_occurrence', 'icd10']), col('BETA_Burden')),
        )
        .withColumn(
            'oddsRatioConfidenceIntervalLower',
            when(
                col('trait_type').isin(['categorical', 'icd_first_occurrence', 'icd10']),
                col('BETA_Burden') - col('SE_Burden'),
            ),
        )
        .withColumn(
            'oddsRatioConfidenceIntervalUpper',
            when(
                col('trait_type').isin(['categorical', 'icd_first_occurrence', 'icd10']),
                col('BETA_Burden') + col('SE_Burden'),
            ),
        )
        .withColumn('studySampleSize', (col('n_cases') + coalesce('n_controls', lit(0))))
        .withColumnRenamed('n_cases', 'studyCases')
        .withColumnRenamed('annotation', 'statisticalMethod')
        .withColumn('statisticalMethodOverview', col('statisticalMethod'))
        .replace(to_replace=METHOD_DESC, subset=['statisticalMethodOverview'])
        .select(to_keep)
        .distinct()
    )


def get_parser():
    'Get parser object for script GenebassGeneBurden.py.'
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--genebass_data',
        help='Input parquet files with Genebass\'s burden associations.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Output gzipped json file following the gene_burden evidence data model.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--log_file',
        help='Destination of the logs generated by this script. Defaults to None',
        type=str,
        nargs='?',
        default=None,
    )

    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()

    # Logger initializer. If no log_file is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    spark = initialize_sparksession()

    evd_df = main(
        genebass_data=args.genebass_data,
        spark_instance=spark,
    )

    write_evidence_strings(evd_df, args.output)
    logging.info(f'Evidence strings have been saved to {args.output}. Exiting.')
