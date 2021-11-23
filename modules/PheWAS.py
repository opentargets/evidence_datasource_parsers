#!/usr/bin/env python

import argparse
import logging
import gzip
import json
import sys

from pyspark import SparkFiles
from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, element_at, split, lit, count, concat, translate
from pyspark.sql.types import IntegerType, DoubleType


class phewasEvidenceGenerator:

    def __init__(self):
        # Create spark session
        sparkConf = (
            SparkConf()
            .set('spark.driver.host', '127.0.0.1')
            .set('spark.driver.memory', '15g')
            .set('spark.executor.memory', '15g')
            .set('spark.driver.maxResultSize', '0')
            .set('spark.debug.maxToStringFields', '2000')
            .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        )
        self.spark = (
            SparkSession.builder
            .config(conf=sparkConf)
            .master('local[*]')
            .getOrCreate()
        )

        # Initialize variables
        self.dataframe = None
        self.enrichedDataframe = None

    def generateEvidenceFromSource(self, inputFile, consequencesFile, diseaseMapping, skipMapping):
        '''
        Processing of the dataframe to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''

        # Read input file
        self.dataframe = (
            self.spark
            .read.csv(inputFile, header=True)
            .select(
                translate(col('gene'), '*', '').alias('gene_symbol'),
                'snp', 'phewas_code', 'phewas_string',
                col('cases').cast(IntegerType()),
                col('odds_ratio').cast(DoubleType()),
                col('p').cast(DoubleType())
            )

            .filter(
                # Filter out nulls and DNA regions. Ex: intergenic, Intergenic, HLA-region
                (col('gene').isNotNull())
                & (~col('gene').contains('ntergenic'))
                & (~col('gene').contains('region'))

                # Keep only significant associations
                & (col('p') < 0.05)
            )
        )

        # Mapping step
        if not skipMapping:
            try:
                self.spark.sparkContext.addFile(diseaseMapping)
                phewasMapping = (
                    self.spark.read.csv(SparkFiles.get(diseaseMapping.split('/')[-1]), sep=r'\t', header=True)
                    .select('Phewas_string', col('EFO_id').alias('EFO_link'))
                    .withColumn('EFO_id', element_at(split(col('EFO_link'), '/'), -1))
                )
                self.dataframe = self.dataframe.join(phewasMapping, on=['Phewas_string'], how='left')
                logging.info('Disease mappings have been imported.')

            except Exception as e:
                logging.error(f'An error occurred while importing disease mappings: \n{e}.')

            else:
                # Filter out invalid disease IDs: MPATH_579, CHEBI_36047
                pattern = r'(^NCIT_C\d+$|^Orphanet_\d+$|^GO_\d+$|^HP_\d+$|^EFO_\d+$|^MONDO_\d+$|^DOID_\d+$|^MP_\d+$)'
                self.dataframe = self.dataframe.filter(col('EFO_id').rlike(pattern))

        else:
            logging.info('Disease mapping has been skipped.')
            self.dataframe = self.dataframe.withColumn('EFO_id', lit(None))

        # Get functional consequence per variant from OT Genetics Portal
        cols = [
            'phewas_string', 'phewas_code', 'EFO_id', 'odds_ratio', 'p',
            'cases', 'gene_symbol', 'consequence_id', 'variantId', 'snp'
        ]
        self.enrichedDataframe = (
            self.enrichVariantData(consequencesFile)
            .dropDuplicates(cols)
        )
        logging.info('Functional consequences have been imported.')

        # Build evidence strings per row
        logging.info('Generating evidence:')
        evidences = (
            self.enrichedDataframe.rdd
            .map(phewasEvidenceGenerator.parseEvidenceString)
            .collect()
        )

        for evidence in evidences:
            if skipMapping:
                # Delete empty keys if mapping is skipped
                del evidence['diseaseFromSourceMappedId']
            if not evidence['variantId']:
                # Delete empty variantId
                del evidence['variantId']

        return evidences

    def enrichVariantData(self, consequencesFile):
        self.spark.sparkContext.addFile(consequencesFile)
        phewasWithConsequences = (
            self.spark.read.csv(SparkFiles.get(consequencesFile.split('/')[-1]), header=True)
            .select(
                col('rsid').alias('snp'),
                col('pos').cast(IntegerType()),
                'gene_symbol', 'chrom', 'ref', 'alt',
                'consequence_link'
            )
            .withColumn(
                'consequence_id',
                element_at(split(col('consequence_link'), '/'), -1)
            )
        )

        # We want to remove all the SNPs associated with many variants
        one2manyVariants = (
            phewasWithConsequences
            .groupBy('snp')
            .agg(count('snp'))
            .filter(col('count(snp)') > 1)
            .toPandas()['snp']
            .tolist()
        )
        phewasWithConsequences = phewasWithConsequences.filter(
            ~col('snp').isin(one2manyVariants)
        )

        # Enriching dataframe with consequences --> more records due to 1:many associations
        self.dataframe = self.dataframe.join(
            phewasWithConsequences,
            on=['gene_symbol', 'snp'],
            how='left'
        )

        # Building variantId: 'chrom_pos_ref_alt' of the respective rsId
        self.dataframe = (self.dataframe.select(
            '*',
            concat(
                col('chrom'),
                lit('_'),
                col('pos'),
                lit('_'),
                col('ref'),
                lit('_'),
                col('alt')
            )
            .alias('variantId')
        ))

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                'datasourceId': 'phewas_catalog',
                'datatypeId': 'genetic_association',
                'diseaseFromSource': row['phewas_string'],
                'diseaseFromSourceId': row['phewas_code'],
                'diseaseFromSourceMappedId': row['EFO_id'],
                'oddsRatio': row['odds_ratio'],
                'resourceScore': row['p'],
                'studyCases': row['cases'],
                'targetFromSourceId': row['gene_symbol'],
                'variantFunctionalConsequenceId': row['consequence_id'] if row['consequence_id'] else 'SO_0001060',
                'variantId': row['variantId'],
                'variantRsId': row['snp']
            }
            return evidence
        except Exception as e:
            raise


def main(inputFile, consequencesFile, diseaseMapping, skipMapping):
    # Initialize evidence builder object
    evidenceBuilder = phewasEvidenceGenerator()

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.generateEvidenceFromSource(inputFile, consequencesFile, diseaseMapping, skipMapping)

    with gzip.open(outputFile, 'wt') as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')

    logging.info(f'{len(evidences)} evidence strings saved into {outputFile}. Exiting.')


if __name__ == '__main__':

    # Initiating parser
    parser = argparse.ArgumentParser(description='This script generates evidences from the PheWAS Catalog data source.')

    parser.add_argument('-i', '--inputFile', required=True, type=str, help='Input .csv file with the table containing association details.')
    parser.add_argument('-c', '--consequencesFile', required=True, type=str, help='Input look-up table containing the variation consequences coming from the Variant Index.')
    parser.add_argument('-d', '--diseaseMapping', required=False, type=str, help='Input look-up table containing the phenotype mappings to an EFO ID.')
    parser.add_argument('-o', '--outputFile', required=True, type=str, help='Name of the compressed json.gz output file containing the evidence strings.')
    parser.add_argument('-s', '--skipMapping', required=False, action='store_true', help='State whether to skip the disease to EFO mapping step.')
    parser.add_argument('-l', '--logFile', help='Destination of the logs generated by this script.', type=str, required=False)

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    consequencesFile = args.consequencesFile
    diseaseMapping = args.diseaseMapping
    outputFile = args.outputFile
    skipMapping = args.skipMapping

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    # Logging parameters
    logging.info(f'PheWAS input table: {inputFile}')
    logging.info(f'Phewas enriched with consequences input file: {consequencesFile}')
    logging.info(f'Phewas phenotype to EFO ID table: {diseaseMapping}')
    logging.info(f'Output file: {outputFile}')

    main(inputFile, consequencesFile, diseaseMapping, skipMapping)
