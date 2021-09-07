#!/usr/bin/env python3
"""Evidence parser for Orphanet's gene-disease associations."""

import argparse
import logging
import os
import sys
import tempfile
import time
from itertools import chain

import xml.etree.ElementTree as ET
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, Row, SparkSession
from pyspark.sql.functions import array_distinct, col, create_map, explode, first, lit, split
from pyspark.sql.types import ArrayType, StringType, StructField, StructType

from ontoma import OnToma

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


class ontoma_efo_lookup():
    """
    Map orphanet diseases to the EFO ontology
    """
    def __init__(self):
        self.otmap = OnToma()

    def get_mapping(self, terms=[]):
        disease_label, disease_id = terms
        label_mapping = [
            result.id_ot_schema
            for result in self.otmap.find_term(disease_label)]

        if len(label_mapping) != 0:
            return (disease_label, disease_id, label_mapping)
        else:
            id_mapping = [
                result.id_ot_schema
                for result in self.otmap.find_term(disease_id)]
            return (disease_label, disease_id, id_mapping)

def parse_orphanet_xml(orphanet_file: str) -> list:
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

    return orphanet_disorders


def main(input_file: str, output_file: str, local: bool = False) -> None:

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

    # Initialize mapping object:
    ol_obj = ontoma_efo_lookup()

    # Map association type to sequence ontology ID:
    so_mapping_expr = create_map([lit(x) for x in chain(*CONSEQUENCE_MAP.items())])

    # Parsing xml file:s
    orphanet_disorders = parse_orphanet_xml(input_file)

    # Create a spark dataframe from the parsed data:
    orphanet_df = (
        spark.createDataFrame(Row(**x) for x in orphanet_disorders)
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
        .persist()
    )

    # Generating a lookup table for the mapped orphanet terms:
    orphanet_diseases = (
        orphanet_df
        .select('diseaseFromSource', 'diseaseFromSourceId')
        .distinct()
        .collect()
    )

    mapped_diseases = [ol_obj.get_mapping(x) for x in orphanet_diseases]
    
    schema = StructType([
        StructField('diseaseFromSource', StringType(), True),
        StructField('diseaseFromSourceId', StringType(), True),
        StructField('diseaseFromSourceMappedId', ArrayType(StringType()), True)
    ])
    mapped_diseases_df = (
        # Dataframe from list of label/id/mapped_id tuples
        spark.createDataFrame(mapped_diseases, schema=schema)
        # Coalesce df row-wise
        .groupBy(
            'diseaseFromSource', 'diseaseFromSourceId')
        .agg(first("diseaseFromSourceMappedId", ignorenulls=True).alias("diseaseFromSourceMappedId"))
        # Explode cases where one diseaseFromSource maps to many EFO terms
        .withColumn('diseaseFromSourceMappedId', explode('diseaseFromSourceMappedId'))
    )

    # Adding EFO mapping as new column and prepare evidence:
    evidence = (
        orphanet_df
        .join(mapped_diseases_df, on=['diseaseFromSource', 'diseaseFromSourceId'], how='left')
        .select(
            'datasourceId', 'datatypeId', 'alleleOrigins', 'confidence', 'diseaseFromSource',
            'diseaseFromSourceId', 'diseaseFromSourceMappedId', 'literature', 'targetFromSource',
            'targetFromSourceId'
        )
    )

    # Save data
    write_evidence_strings(evidence, output_file)
    logging.info(f'{orphanet_df.count()} evidence strings have been saved to {output_file}')

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
    parser.add_argument(
        '--local', action='store_true', required=False, default=False,
        help='Flag to indicate if the script is executed locally or on the cluster'
    )

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    log_file = args.logFile
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

    main(input_file, output_file, is_local)
