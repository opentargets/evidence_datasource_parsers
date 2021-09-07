#!/usr/bin/env python3
"""Evidence parser for ClinGen's Gene Validity Curations."""

import argparse
import json
import logging
import os
import tempfile

from ontoma import OnToma
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import array, col, explode, first, lit, monotonically_increasing_id, struct, trim
from pyspark.sql.types import ArrayType, StringType, StructField, StructType, TimestampType


class ontoma_efo_lookup():
    """
    Map orphanet diseases to the EFO ontology
    """
    def __init__(self):
        self.otmap = OnToma(cache_dir='cache')

    def get_mapping(self, terms=[]):
        disease_label, disease_id = terms
        label_mapping = [
            result.id_ot_schema
            for result in self.otmap.find_term(disease_label)]

        if len(label_mapping) != 0:
            return (disease_label, disease_id, label_mapping)
        try:
            id_mapping = [
                result.id_ot_schema
                for result in self.otmap.find_term(disease_id)]
            return (disease_label, disease_id, id_mapping)
        except AttributeError:
            return None

def main(input_file: str, output_file: str, local: bool = False):

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

    # Read Gene Validity Curations input file
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
        spark.read.csv(input_file, schema=clingen_schema)
        # The first 6 rows of the GVC file include metadata that needs to be dropped
        .withColumn('idx', monotonically_increasing_id())
        .filter(col('idx') > 5)
        .drop('idx')
    )

    # Generating a lookup table for the mapped orphanet terms:
    diseases = (
        clingen_df
        .select('DISEASE LABEL', 'DISEASE ID (MONDO)')
        .distinct()
        .collect()
    )

    mapped_diseases = [ol_obj.get_mapping(x) for x in diseases]

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

    # Adding EFO mapping as new column:
    clingen_df = (
        clingen_df
        .withColumnRenamed('DISEASE LABEL', 'diseaseFromSource')
        .withColumnRenamed('DISEASE ID (MONDO)', 'diseaseFromSourceId')
        .join(mapped_diseases_df, on=['diseaseFromSource', 'diseaseFromSourceId'], how='left')
    )

    # Turn table into evidence strings
    evidence = process_clingen(clingen_df)
    write_evidence_strings(evidence, output_file)
    logging.info(f'{evidence.count()} evidence strings have been saved to {output_file}')


def process_clingen(clingen_df: DataFrame):
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
            # col('DISEASE LABEL').alias('diseaseFromSource'),
            # col('DISEASE ID (MONDO)').alias('diseaseFromSourceId'),
            array(col('MOI')).alias('allelicRequirements'),
            array(col('urls')).alias('urls'),
            col('CLASSIFICATION').alias('confidence'),
            col('GCEP').alias('studyId'),
            'diseaseFromSource', 'diseaseFromSourceId',
            'diseaseFromSourceMappedId'
        )
    )

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


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse ClinGen gene-disease associations from Gene Validity Curations')
    parser.add_argument('-i', '--input_file',
                        help='Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Absolute path of the gzipped, JSON evidence file.',
                        type=str, required=True)
    parser.add_argument('-l', '--log_file', type=str,
                        help='Optional filename to redirect the logs into.')

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

    main(args.input_file, args.output_file)
