#!/usr/bin/env python3
"""Evidence parser for Orphanet's gene-disease associations."""

import argparse
from itertools import chain
import logging
import os
import sys
import tempfile

import xml.etree.ElementTree as ET
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, Row, SparkSession
from pyspark.sql.functions import array_distinct, col, create_map, lit, split

from common.ontology import add_efo_mapping

# The rest of the types are assigned to -> germline for allele origins
EXCLUDED_ASSOCIATIONTYPES = [
    "Major susceptibility factor in",
    "Part of a fusion gene in",
    "Candidate gene tested in",
    "Role in the phenotype of",
    "Biomarker tested in",
    "Disease-causing somatic mutation(s) in"
]

# Assigning variantFunctionalConsequenceId:
CONSEQUENCE_MAP = {
    "Disease-causing germline mutation(s) (loss of function) in": "SO_0002054",
    "Disease-causing germline mutation(s) in": None,
    "Modifying germline mutation in": None,
    "Disease-causing germline mutation(s) (gain of function) in": "SO_0002053"
}


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

    # Read and process Orphanet's XML file into evidence strings
    
    orphanet_df = parse_orphanet_xml(input_file, spark)
    logging.info('Orphanet input file has been imported. Processing evidence strings.')

    evidence_df = process_orphanet(orphanet_df)

    evidence_df = add_efo_mapping(evidence_strings=evidence_df, spark_instance=spark, ontoma_cache_dir=cache_dir)
    logging.info('Disease mappings have been added.')

    # Save data
    write_evidence_strings(evidence_df, output_file)
    logging.info(f'{evidence_df.count()} evidence strings have been saved to {output_file}')

def parse_orphanet_xml(orphanet_file: str, spark_instance) -> DataFrame:
    """
    Function to parse Orphanet xml dump and return the parsed
    data as a list of dictionaries.

    Args:
        orphanet_file (str): Orphanet XML filename

    Returns:
        parsed data as a list of dictionaries
    """

    # Reading + validating xml:
    tree = ET.parse(orphanet_file)
    assert isinstance(tree, ET.ElementTree)

    root = tree.getroot()
    assert isinstance(root, ET.Element)

    # Checking if the basic nodes are in the xml structure:

    logging.info(f"There are {root.find('DisorderList').get('count')} disease in the Orphanet xml file.")

    orphanet_disorders = []

    for disorder in root.find('DisorderList').findall('Disorder'):

        # Extracting disease information:
        parsed_disorder = {
            "diseaseFromSource": disorder.find('Name').text,
            "diseaseFromSourceId": 'Orphanet_' + disorder.find('OrphaCode').text,
            "type": disorder.find('DisorderType/Name').text,
        }

        # One disease might be mapped to multiple genes:
        for association in disorder.find('DisorderGeneAssociationList'):

            # For each mapped gene, an evidence is created:
            evidence = parsed_disorder.copy()

            # Not all gene/disease association is backed up by publication:
            try:
                evidence['literature'] = [
                    pmid.replace('[PMID]', '').rstrip() for pmid in association.find('SourceOfValidation').text.split('_') if '[PMID]' in pmid
                ]
            except AttributeError:
                evidence['literature'] = None

            evidence['associationType'] = association.find('DisorderGeneAssociationType/Name').text
            evidence['confidence'] = association.find('DisorderGeneAssociationStatus/Name').text

            # Parse gene name and id - going for Ensembl gene id only:
            gene = association.find('Gene')
            evidence['targetFromSource'] = gene.find('Name').text

            # Extracting ensembl gene id from cross references:
            ensembl_gene_id = [
                xref.find('Reference').text for xref in gene.find('ExternalReferenceList') if 'ENSG' in xref.find('Reference').text
            ]
            evidence['targetFromSourceId'] = ensembl_gene_id[0] if len(ensembl_gene_id) > 0 else None

            # Collect evidence:
            orphanet_disorders.append(evidence)

    # Create a spark dataframe from the parsed data
    orphanet_df = spark_instance.createDataFrame(Row(**x) for x in orphanet_disorders)

    return orphanet_df

def process_orphanet(orphanet_df: DataFrame) -> DataFrame:
    """
    The JSON Schema format is applied to the df
    """

    # Map association type to sequence ontology ID:
    so_mapping_expr = create_map([lit(x) for x in chain(*CONSEQUENCE_MAP.items())])

    evidence_df = (
        orphanet_df
        .filter(
            ~col('associationType').isin(EXCLUDED_ASSOCIATIONTYPES)
        )
        .filter(~col('targetFromSourceId').isNull())
        .withColumn('dataSourceId', lit('orphanet'))
        .withColumn('datatypeId', lit('genetic_association'))
        .withColumn('alleleOrigins', split(lit('germline'), "_"))
        .withColumn('literature', array_distinct(col('literature')))
        .withColumn('variantFunctionalConsequenceId', so_mapping_expr.getItem(col('associationType')))
        .drop('associationType', 'type')

        # Select the evidence relevant fields
        .select(
            'datasourceId', 'datatypeId', 'alleleOrigins', 'confidence', 'diseaseFromSource',
            'diseaseFromSourceId', 'literature', 'targetFromSource',
            'targetFromSourceId'
        )
        .persist()
    )

    return evidence_df


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=('Parse Orphanet gene-disease annotation downloaded from '
                                                  'http://www.orphadata.org/data/xml/en_product6.xml'))
    parser.add_argument('--input_file', help='Xml file containing target/disease associations.', type=str)
    parser.add_argument('--output_file', help='Absolute path of the gzipped, JSON evidence file.', type=str)
    parser.add_argument('--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    parser.add_argument('--cache_dir', required=False, help='Directory to store the OnToma cache files in.')
    parser.add_argument(
        '--local', action='store_true', required=False, default=False,
        help='Flag to indicate if the script is executed locally or on the cluster'
    )

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    log_file = args.logFile
    cache_dir = args.cache_dir
    is_local = args.local

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    main(input_file, output_file, cache_dir, is_local)
