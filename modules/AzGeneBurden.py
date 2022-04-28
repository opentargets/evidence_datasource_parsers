#!/usr/bin/env python
"""This module extracts and processes target/disease evidence from the AstraZeneca PheWAS Portal."""

import argparse
import logging
import sys

from pandas import read_excel
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame
from pyspark.sql.functions import array, col, explode, expr, lit, log10, pow, round, when
from pyspark.sql.types import DoubleType, IntegerType

from common.evidence import initialize_sparksession, write_evidence_strings

METHOD_DESC = {
    'ptv': 'Burden test carried out with PTVs with a MAF smaller than 0.1%.',
    'ptv5pcnt': 'Burden test carried out with PTVs with a MAF smaller than 5%.',
    'UR': 'Burden test carried out with ultra rare damaging variants (MAF ≈ 0%).',
    'URmtr': 'Burden test carried out with MTR-informed ultra rare damaging variants (MAF ≈ 0%).',
    'raredmg': 'Burden test carried out with rare missense variants with a MAF smaller than 0.005%.',
    'raredmgmtr': 'Burden test carried out with MTR-informed rare missense variants with a MAF smaller than 0.005%.',
    'flexdmg': 'Burden test carried out with damaging variants with a MAF smaller than 0.01%.',
    'flexnonsyn': 'Burden test carried out with non synonymous variants with a MAF smaller than 0.01%.',
    'flexnonsynmtr': 'Burden test carried out with MTR-informed non synonymous variants with a MAF smaller than 0.01%.',
    'ptvraredmg': 'Burden test carried out with PTV or rare missense variants.',
    'rec': 'Burden test carried out with non-synonymous recessive variants with a MAF smaller than 1%.',
}


def main(az_binary_data: str, az_quant_data: str, az_trait_mappings: str, spark_instance: SparkSession) -> DataFrame:
    """
    This module extracts and processes target/disease evidence from the raw AstraZeneca PheWAS Portal.
    """
    logging.info(f'File with the AZ PheWAS Portal binary traits associations: {az_binary_data}')
    logging.info(f'File with the AZ PheWAS Portal quantitative traits associations: {az_quant_data}')

    # Load data
    az_trait_mappings_df = (
        spark_instance.createDataFrame(
            read_excel(az_trait_mappings, sheet_name='Sheet1').filter(items=['Phenotype', 'EFO'])
        )
        .withColumn('EFO', explode(expr("regexp_extract_all(EFO, '([A-Za-z]+_[0-9]+)')")))
        .distinct()
    )
    az_phewas_df = (
        spark_instance.read.parquet(az_binary_data)
        # Renaming of some columns to match schemas of both binary and quantitative evidence
        .withColumnRenamed('BinOddsRatioLCI', 'LCI')
        .withColumnRenamed('BinOddsRatioUCI', 'UCI')
        .withColumnRenamed('BinNcases', 'nCases')
        .withColumnRenamed('BinQVcases', 'nCasesQV')
        .withColumnRenamed('BinNcontrols', 'nControls')
        # Combine binary and quantitative evidence into one dataframe
        .unionByName(
            spark_instance.read.parquet(az_quant_data)
            .withColumn('nCases', col('nSamples'))
            .withColumnRenamed('YesQV', 'nCasesQV'),
            allowMissingColumns=True,
        )
        # Add mapped trait from GWAS Catalog WIP sheet
        # This is a temporary solution, as the AZ studies are not yet integrated into the GWAS Catalog
        .join(az_trait_mappings_df, on='Phenotype', how='left')
        .withColumn('pValue', col('pValue').cast(DoubleType()))
        .filter(col('pValue') <= 1e-7)
        .distinct()
        .repartition(20)
        .persist()
    )

    # Filter out associations from the synonymous model used as a negative control
    az_phewas_df = remove_false_positives(az_phewas_df)

    # WARNING: There are some associations with a p-value of 0.0 in the AstraZeneca PheWAS Portal.
    # This is a bug we still have to ellucidate and it might be due to a float overflow.
    # These evidence need to be manually corrected in order not to lose them and for them to pass validation
    # As an interim solution, their p value will equal to the minimum in the evidence set.
    logging.warning(f"There are {az_phewas_df.filter(col('pValue') == 0.0).count()} evidence with a p-value of 0.0.")
    minimum_pvalue = az_phewas_df.filter(col('pValue') > 0.0).agg({'pValue': 'min'}).collect()[0]['min(pValue)']
    az_phewas_df = az_phewas_df.withColumn(
        'pValue', when(col('pValue') == 0.0, minimum_pvalue).otherwise(col('pValue'))
    )

    # Write output
    evd_df = parse_az_phewas_evidence(az_phewas_df)

    if evd_df.filter(col('resourceScore') == 0).count() != 0:
        logging.error('There are evidence with a P value of 0.')
        raise AssertionError(
            f"There are {evd_df.filter(col('resourceScore') == 0).count()} evidence with a P value of 0."
        )

    if not 17000 < evd_df.count() < 18000:
        logging.error(f'AZ PheWAS Portal number of evidence are different from expected: {evd_df.count()}')
        raise AssertionError('AZ PheWAS Portal number of evidence are different from expected.')
    logging.info(f'{evd_df.count()} evidence strings have been processed.')

    return evd_df


def remove_false_positives(az_phewas_df: DataFrame) -> DataFrame:
    """Remove associations present in the synonymous negative control."""

    false_positives = az_phewas_df.filter(col('CollapsingModel') == 'syn').select('Gene', 'Phenotype').distinct()
    true_positives = az_phewas_df.join(false_positives, on=['Gene', 'Phenotype'], how='left_anti').distinct()
    logging.info(
        f'{az_phewas_df.count() - true_positives.count()} false positive evidence of association have been dropped.'
    )

    return true_positives


def parse_az_phewas_evidence(az_phewas_df: DataFrame) -> DataFrame:
    """
    Parse Astra Zeneca's PheWAS Portal evidence.
    Args:
        az_phewas_df: DataFrame with Astra Zeneca's PheWAS Portal data
    Returns:
        evd_df: DataFrame with Astra Zeneca's data following the t/d evidence schema.
    """
    to_keep = [
        'datasourceId',
        'datatypeId',
        'targetFromSourceId',
        'diseaseFromSource',
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
        'literature',
        'projectId',
        'cohortId',
        'studySampleSize',
        'studyCases',
        'studyCasesWithQualifyingVariants',
        'statisticalMethod',
        'statisticalMethodOverview',
    ]

    return (
        az_phewas_df.withColumn('datasourceId', lit('gene_burden'))
        .withColumn('datatypeId', lit('genetic_association'))
        .withColumn('literature', array(lit('34375979')))
        .withColumn('projectId', lit('AstraZeneca PheWAS Portal'))
        .withColumn('cohortId', lit('UK Biobank 450k'))
        .withColumnRenamed('Gene', 'targetFromSourceId')
        .withColumnRenamed('Phenotype', 'diseaseFromSource')
        .withColumnRenamed('EFO', 'diseaseFromSourceMappedId')
        .withColumn('resourceScore', col('pValue'))
        .withColumn('pValueExponent', log10(col('pValue')).cast(IntegerType()) - lit(1))
        .withColumn('pValueMantissa', round(col('pValue') / pow(lit(10), col('pValueExponent')), 3))
        .withColumn(
            'beta',
            when(col('Type') == 'Quantitative', col('beta')),
        )
        .withColumn(
            'betaConfidenceIntervalLower',
            when(col('Type') == 'Quantitative', col('LCI')),
        )
        .withColumn(
            'betaConfidenceIntervalUpper',
            when(col('Type') == 'Quantitative', col('UCI')),
        )
        .withColumn(
            'oddsRatio',
            when(col('Type') == 'Binary', col('binOddsRatio')),
        )
        .withColumn(
            'oddsRatioConfidenceIntervalLower',
            when(col('Type') == 'Binary', col('LCI')),
        )
        .withColumn(
            'oddsRatioConfidenceIntervalUpper',
            when(col('Type') == 'Binary', col('UCI')),
        )
        .withColumn('ancestry', lit('EUR'))
        .withColumn('ancestryId', lit('HANCESTRO_0005'))
        .withColumnRenamed('nSamples', 'studySampleSize')
        .withColumnRenamed('nCases', 'studyCases')
        .withColumnRenamed('nCasesQV', 'studyCasesWithQualifyingVariants')
        .withColumnRenamed('CollapsingModel', 'statisticalMethod')
        .withColumn('statisticalMethodOverview', col('statisticalMethod'))
        .replace(to_replace=METHOD_DESC, subset=['statisticalMethodOverview'])
        .select(to_keep)
        .distinct()
    )


def get_parser():
    'Get parser object for script AzGeneBurden.py.'
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--az_binary_data',
        help='Input parquet files with AZ\'s PheWAS associations of binary traits.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--az_quant_data',
        help='Input parquet files with AZ\'s PheWAS associations of quantitative traits.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--az_trait_mappings',
        help='Input Excel containing the AZ traits with their EFO mappings.',
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
        az_binary_data=args.az_binary_data,
        az_quant_data=args.az_quant_data,
        az_trait_mappings=args.az_trait_mappings,
        spark_instance=spark,
    )

    write_evidence_strings(evd_df, args.output)
    logging.info(f'Evidence strings have been saved to {args.output}. Exiting.')
