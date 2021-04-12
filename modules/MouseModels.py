#!/usr/bin/env python3

import argparse
import collections
import datetime
import gzip
import hashlib
import json
import logging
import os
import pathlib
import re
import urllib.request

import pyspark
import pyspark.sql.functions
import requests
from retry import retry

from common.RareDiseasesUtils import RareDiseaseMapper
from ontologyutils.rdf_utils import OntologyClassReader
import opentargets.model.core as cttv
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.phenotype as evidence_phenotype
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.association_score as association_score
from settings import Config

HGNC_DATASET_URI = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
HGNC_DATASET_FILENAME = os.path.split(HGNC_DATASET_URI)[-1]
MGI_DATASET_URI = 'http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt'
MGI_DATASET_FILENAME = os.path.split(MGI_DATASET_URI)[-1]

IMPC_SOLR_HOST = 'http://www.ebi.ac.uk/mi/impc/solr/phenodigm/select'
IMPC_SOLR_DATA_TYPES = ('gene', 'gene_gene', 'mouse_model', 'disease_model_summary', 'disease_gene_summary', 'disease',
                        'ontology_ontology', 'ontology')
IMPC_FILENAME = 'impc_solr_{data_type}.json'


class ImpcSolrRetriever:
    """Retrieves data from the IMPC SOLR API and saves JSONs to the specified location."""

    def __init__(self, solr_host, timeout, rows):
        """Initialise query parameters.

        Args:
            solr_host (str): SOLR endpoint to make the requests against.
            timeout (int): Timeout to apply to the requests, in seconds.
            rows (int): Number of SOLR documents returned in a single batch."""
        self.solr_host = solr_host
        self.timeout = timeout
        self.rows = rows

    # The @retry decorator ensures that the requests are retried in case of network or server error.
    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def query_solr(self, data_type, start):
        """Returns one batch of the SOLR documents of the specified data type."""
        params = {'q': '*:*', 'fq': f'type:{data_type}', 'start': start, 'rows': self.rows}
        response = requests.get(self.solr_host, params=params, timeout=self.timeout)
        # Check for HTTP errors
        response.raise_for_status()
        return response.json()

    def fetch_data(self, data_type, output_filename):
        """Fetch all data of the requested type to the specified location."""
        output_file = open(output_filename, 'wt')
        # Initialise the counters.
        start, total = 0, 0
        while True:
            solr_data = self.query_solr(data_type, start)
            # If we don't find any items, we break.
            if solr_data['response']['numFound'] == 0:
                break
            # Write data to file.
            for doc in solr_data['response']['docs']:
                json.dump(doc, output_file)
                output_file.write('\n')
            # Increment the counters.
            start += self.rows
            total += len(solr_data['response']['docs'])
            # Exit if all documents have been retrieved.
            if total == solr_data['response']['numFound']:
                break
        output_file.close()


class PhenoDigm(RareDiseaseMapper):

    def __init__(self, logging):
        super(PhenoDigm, self).__init__()
        self.spark = pyspark.sql.SparkSession.builder.appName('phenodigm_parser').getOrCreate()
        self.hgnc_id_to_ensembl_id, self.mgi_id_to_ensembl_id = None, None
        self.disease_model_summary = None
        self.evidence = None
        self.efo = OntologyClassReader()
        self.hpo = OntologyClassReader()
        self.mp = OntologyClassReader()
        self.cache = collections.OrderedDict()
        self.counter = 0
        self._logger = logging

    def update_cache(self, cache_dir):
        """Fetch the Ensembl gene ID and SOLR data into the local cache directory."""
        pathlib.Path(cache_dir).mkdir(parents=False, exist_ok=True)

        self._logger.info('Fetching human gene ID mappings from HGNC.')
        urllib.request.urlretrieve(HGNC_DATASET_URI, os.path.join(cache_dir, HGNC_DATASET_FILENAME))

        self._logger.info('Fetching mouse gene ID mappings from MGI.')
        urllib.request.urlretrieve(MGI_DATASET_URI, os.path.join(cache_dir, MGI_DATASET_FILENAME))

        self._logger.info('Fetching Phenodigm data from IMPC SOLR.')
        impc_solr_retriever = ImpcSolrRetriever(solr_host=IMPC_SOLR_HOST, timeout=600, rows=50000)
        for data_type in IMPC_SOLR_DATA_TYPES:
            self._logger.info(f'Fetching Phenodigm data type {data_type}.')
            filename = os.path.join(cache_dir, IMPC_FILENAME.format(data_type=data_type))
            impc_solr_retriever.fetch_data(data_type, filename)

    def load_data_from_cache(self, cache_dir):
        """Load the Ensembl gene ID and SOLR data from the downloaded files into Spark."""
        self.hgnc_id_to_ensembl_id = (
            self.spark.read.csv(os.path.join(cache_dir, HGNC_DATASET_FILENAME), sep='\t', header=True)
                .select('hgnc_id', 'ensembl_gene_id')
            .withColumnRenamed('ensembl_gene_id', 'ensembl_human_gene_id')
        )
        self.mgi_id_to_ensembl_id = (
            self.spark.read.csv(os.path.join(cache_dir, MGI_DATASET_FILENAME), sep='\t', header=True)
                .withColumnRenamed('1. MGI accession id', 'mgi_id')
                .withColumnRenamed('11. Ensembl gene id', 'ensembl_mouse_gene_id')
                .select('mgi_id', 'ensembl_mouse_gene_id')
        )
        self.disease_model_summary = (
            self.spark.read.json(os.path.join(cache_dir, IMPC_FILENAME.format(data_type='disease_model_summary')))
                .select('model_id', 'model_genetic_background', 'model_description', 'disease_id', 'disease_term',
                        'disease_model_max_norm')
                .limit(100)
        )

    def generate_phenodigm_evidence_strings(self):
        self.evidence = (
            self.disease_model_summary

            # Rename some columns
            .withColumnRenamed('disease_id', 'diseaseFromSourceId')
            .withColumnRenamed('disease_term', 'diseaseFromSource')
            .withColumnRenamed('model_description', 'biologicalModelAllelicComposition')
            .withColumnRenamed('model_genetic_background', 'biologicalModelGeneticBackground')

            # Strip trailing modifiers from the model ID
            # For example: MGI:6274930#hom#early → MGI:6274930
            .withColumn(
                'biologicalModelId',
                pyspark.sql.functions.split(pyspark.sql.functions.col('model_id'), '#').getItem(0)
            )

            # Convert the percentage score into fraction
            .withColumn('resourceScore', pyspark.sql.functions.col('disease_model_max_norm') / 100.0)

            # Remove intermediate columns
            .drop('disease_model_max_norm', 'model_id')

            # Add constant value columns
            .withColumn('datasourceId', pyspark.sql.functions.lit('phenodigm'))
            .withColumn('datatypeId', pyspark.sql.functions.lit('animal_model'))
        )

    def write_evidence_strings(self, filename):
        self.evidence.write.format('json').save(filename)

    def process_ontologies(self):
        self._logger.info('Load MP classes.')
        self.mp.load_mp_classes()
        self._logger.info('Load HP classes.')
        self.hpo.load_hpo_classes()
        self._logger.info('Load EFO classes.')
        self.efo.load_efo_classes()

    def get_ontology_mappings(self):
        self._logger.info('Get ontology mappings.')
        self.get_omim_to_efo_mappings()
        self._logger.info('Get all ZOOMA mappings for Open Targets.')
        self.get_opentargets_zooma_to_efo_mappings()

    def process_all(self, cache_dir, output, use_cached):
        if not use_cached:
            self._logger.info('Update the gene mapping and SOLR cache')
            self.update_cache(cache_dir)

        self._logger.info('Load gene mappings and SOLR data from local cache')
        self.load_data_from_cache(cache_dir)

        self._logger.info('Process ontologies.')
        # TODO bring back
        # self.process_ontologies()
        # self.get_ontology_mappings()

        self._logger.info('Build evidence strings.')
        self.generate_phenodigm_evidence_strings()

        self._logger.info('Write evidence strings.')
        self.write_evidence_strings(output)

        return


def main(cache_dir, output, use_cached=False, log_file=None):
    # Initialize logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    # Process the data.
    PhenoDigm(logging).process_all(cache_dir, output, use_cached)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evidence parser for animal models sources from PhenoDigm.')
    parser.add_argument('--cache-dir', help='Directory to store intermediate cache files in.', required=True)
    parser.add_argument('--output', help='Name of the JSON file to output the evidence strings.', required=True)
    parser.add_argument('--use-cached', help='Use existing downloaded cache.', action='store_true')
    parser.add_argument('--log-file', help='Optional filename to redirect the logs into.')
    args = parser.parse_args()
    main(args.cache_dir, args.output, args.use_cached, args.log_file)
