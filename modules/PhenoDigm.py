#!/usr/bin/env python3
"""Evidence parser for the animal model sources from PhenoDigm."""

import argparse
import json
import logging
import os
import pathlib
import shutil
import tempfile
import urllib.request

import pyspark
import pyspark.sql.functions as pf
import requests
from retry import retry


# Human gene mappings.
HGNC_DATASET_URI = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
HGNC_DATASET_FILENAME = os.path.split(HGNC_DATASET_URI)[-1]

# Mouse gene mappings.
MGI_DATASET_URI = 'http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt'
MGI_DATASET_FILENAME = os.path.split(MGI_DATASET_URI)[-1]

# Mouse model data from IMPC SOLR.
IMPC_SOLR_HOST = 'http://www.ebi.ac.uk/mi/impc/solr/phenodigm/select'
# Other tables (not currently used): gene, disease_gene_summary.
IMPC_SOLR_TABLES = ('gene_gene', 'mouse_model', 'disease_model_summary', 'disease', 'ontology', 'ontology_ontology')
IMPC_FILENAME = 'impc_solr_{data_type}.json'
IMPC_SOLR_BATCH_SIZE = 100000
IMPC_SOLR_TIMEOUT = 600
DEFAULT_ASSOCIATION_SCORE_CUTOFF = 90.0


class ImpcSolrRetriever:
    """Retrieve data from the IMPC SOLR API and save the JSONs to the specified location."""

    def __init__(self, solr_host: str, timeout: int, rows: int):
        """Initialise the query parameters: SOLR endpoint to make the requests against; timeout to apply to the
        requests, in seconds; and the number of SOLR documents requested in a single batch."""
        self.solr_host = solr_host
        self.timeout = timeout
        self.rows = rows

    # The decorator ensures that the requests are retried in case of network or server errors.
    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def query_solr(self, data_type, start):
        """Request one batch of SOLR documents of the specified data type."""
        params = {'q': '*:*', 'fq': f'type:{data_type}', 'start': start, 'rows': self.rows}
        response = requests.get(self.solr_host, params=params, timeout=self.timeout)
        response.raise_for_status()  # Check for HTTP errors. This will be caught by @retry.
        return response.json()

    def fetch_data(self, data_type, output_filename):
        """Fetch all rows of the requested data type to the specified location."""
        with open(output_filename, 'wt') as outfile:
            start, total = 0, 0  # Initialise the counters.
            while True:
                solr_data = self.query_solr(data_type, start)
                assert solr_data['response']['numFound'] != 0, f'SOLR did not return any data for {data_type}.'
                for doc in solr_data['response']['docs']:  # Write data to file.
                    json.dump(doc, outfile)
                    outfile.write('\n')
                # Increment the counters.
                start += self.rows
                total += len(solr_data['response']['docs'])
                # Exit when all documents have been retrieved.
                if total == solr_data['response']['numFound']:
                    break


class PhenoDigm:
    """Retrieve the data, load it into Spark, process and write the resulting evidence strings."""

    def __init__(self, logger, cache_dir):
        super(PhenoDigm, self).__init__()
        self.logger = logger
        self.cache_dir = cache_dir
        self.spark = pyspark.sql.SparkSession.builder.appName('phenodigm_parser').getOrCreate()
        self.hgnc_gene_id_to_ensembl_human_gene_id, self.mgi_gene_id_to_ensembl_mouse_gene_id = [None] * 2
        self.mouse_gene_to_human_gene, self.mouse_phenotype_to_human_phenotype = [None] * 2
        self.disease_model_summary, self.mouse_model, self.disease, self.ontology = [None] * 4
        self.evidence = None

    def update_cache(self):
        """Fetch the Ensembl gene ID and SOLR data into the local cache directory."""
        pathlib.Path(self.cache_dir).mkdir(parents=False, exist_ok=True)

        self.logger.info('Fetching human gene ID mappings from HGNC.')
        urllib.request.urlretrieve(HGNC_DATASET_URI, os.path.join(self.cache_dir, HGNC_DATASET_FILENAME))

        self.logger.info('Fetching mouse gene ID mappings from MGI.')
        urllib.request.urlretrieve(MGI_DATASET_URI, os.path.join(self.cache_dir, MGI_DATASET_FILENAME))

        self.logger.info('Fetching PhenoDigm data from IMPC SOLR.')
        impc_solr_retriever = ImpcSolrRetriever(solr_host=IMPC_SOLR_HOST, timeout=IMPC_SOLR_TIMEOUT,
                                                rows=IMPC_SOLR_BATCH_SIZE)
        for data_type in IMPC_SOLR_TABLES:
            self.logger.info(f'Fetching PhenoDigm data type {data_type}.')
            filename = os.path.join(self.cache_dir, IMPC_FILENAME.format(data_type=data_type))
            impc_solr_retriever.fetch_data(data_type, filename)

    def load_tsv(self, filename):
        return self.spark.read.csv(os.path.join(self.cache_dir, filename), sep='\t', header=True)

    def load_solr_json(self, data_type):
        return self.spark.read.json(os.path.join(self.cache_dir, IMPC_FILENAME.format(data_type=data_type)))

    def load_data_from_cache(self):
        """Load the Ensembl gene ID and SOLR data from the downloaded TSV/JSON files into Spark."""
        # Mappings from HGNC/MGI gene IDs to Ensembl gene IDs.
        self.hgnc_gene_id_to_ensembl_human_gene_id = (  # E.g. 'HGNC:5', 'ENSG00000121410'.
            self.load_tsv(HGNC_DATASET_FILENAME)
            .withColumnRenamed('hgnc_id', 'hgnc_gene_id')
            .withColumnRenamed('ensembl_gene_id', 'targetFromSourceId')  # Using the final name.
            .select('hgnc_id', 'targetFromSourceId')
        )
        self.mgi_gene_id_to_ensembl_mouse_gene_id = (  # E.g. 'MGI:87853', 'ENSMUSG00000027596'.
            self.load_tsv(MGI_DATASET_FILENAME)
            .withColumnRenamed('1. MGI accession id', 'mgi_gene_id')
            .withColumnRenamed('11. Ensembl gene id', 'targetInModel')  # Using the final name.
            .select('mgi_gene_id', 'targetInModel')
        )

        # Mouse to human mappings.
        self.mouse_gene_to_human_gene = (  # E.g. 'MGI:1346074', 'HGNC:4024'.
            self.load_solr_json('gene_gene')
            .withColumnRenamed('gene_id', 'mgi_gene_id')
            .select('mgi_gene_id', 'hgnc_gene_id')
        )
        self.mouse_phenotype_to_human_phenotype = (  # E. g. 'MP:0000745','HP:0100033'
            self.load_solr_json('ontology_ontology')
            .select('mp_id', 'hp_id')
        )

        # Mouse model and disease data.
        # Note that the models are accessioned with the same prefix ('MGI:') as genes, but they are separate entities.
        self.mouse_model = (  # E. g. 'MGI:3800884', ['MP:0001304 cataract'].
            self.load_solr_json('mouse_model')
            .select('model_id', 'model_phenotypes')
        )
        self.disease = (  # E.g. 'OMIM:609258', ['HP:0000545 Myopia'].
            self.load_solr_json('disease')
            .select('disease_id', 'disease_phenotypes')
        )
        self.disease_model_summary = (
            # E. g. 'MGI:2681494', 'C57BL/6JY-smk', 'smk/smk', 'ORPHA:3097', 'Meacham Syndrome', 91.6, 'MGI:98324'.
            self.load_solr_json('disease_model_summary')
            .withColumnRenamed('marker_id', 'mgi_gene_id')
            .withColumnRenamed('model_description', 'biologicalModelAllelicComposition')  # Using the final name.
            .withColumnRenamed('model_genetic_background', 'biologicalModelGeneticBackground')  # Using the final name.
            .withColumnRenamed('disease_model_max_norm', 'resourceScore')  # Using the final name.
            .select('model_id', 'biologicalModelGeneticBackground', 'biologicalModelAllelicComposition', 'disease_id',
                    'disease_term', 'resourceScore', 'mgi_gene_id')
        )
        self.ontology = (
            self.load_solr_json('ontology')
            .filter((pf.col('ontology') == 'MP') | (pf.col('ontology') == 'HP'))
            .select('phenotype_id', 'phenotype_term')
            .filter(pf.col('on'))
        )
        assert self.ontology.select('phenotype_id').distinct().count() == self.ontology.count(), \
            f'Encountered multiple names for the same term in the ontology table.'

    def generate_phenodigm_evidence_strings(self, score_cutoff):
        """Generate the evidence by renaming, transforming and joining the columns."""
        # Process ontology information to enable MP and HP term lookup based on the ID.
        mp_terms, hp_terms = (
            self.ontology
            .filter(self.ontology == ontology_name)
            .withColumn('phenotype_id', f'{ontology_name}_id')
            .withColumn('phenotype_term', f'{ontology_name}_term')
            for ontology_name in ('mp', 'hp')
        )

        # Split lists of phenotypes in the `mouse_model` and `disease` tables and keep only the ID. For example, one row
        # with ['MP:0001529 abnormal vocalization', 'MP:0002981 increased liver weight'] becomes two rows with
        # 'MP:0001529' and 'MP:0002981'.
        model_phenotypes_split = (
            self.mouse_model
            .withColumn('phenotype', pf.explode('model_phenotypes'))
            .withColumn('mp_id', pf.split(pf.col('phenotype'), ' ').getItem(0))
            .select('model_id', 'mp_id')
        )
        human_phenotypes_split = (
            self.disease
            .withColumn('phenotype', pf.explode('disease_phenotypes'))
            .withColumn('hp_id', pf.split(pf.col('phenotype'), ' ').getItem(0))
            .select('disease_id', 'hp_id')
        )

        # We are reporting all mouse phenotypes for a model, regardless of whether they can be mapped into any human
        # disease.
        all_mouse_phenotypes = (
            model_phenotypes_split
            .join(mp_terms, on='mp_id', how='inner')
            .groupby('model_id')
            .agg(
                pf.collect_set(pf.struct(
                    pf.col('mp_id').alias('id'),
                    pf.col('mp_term').alias('label')
                )).alias('diseaseModelAssociatedModelPhenotypes')
            )
            .select('model_id', 'diseaseModelAssociatedModelPhenotypes')
        )
        # For human phenotypes, we only want to include the ones which are present in the disease *and* also can be
        # traced back to the model phenotypes through the MP → HP mapping relationship.
        matched_human_phenotypes = (
            model_phenotypes_split
            .join(self.mouse_phenotype_to_human_phenotype, on='mp_id', how='inner')
            .join(human_phenotypes_split, on='hp_id', how='inner')
            .join(hp_terms, on='hp_id', how='inner')
            .groupby('model_id', 'disease_id')
            .agg(
                pf.collect_set(pf.struct(
                    pf.col('hp_id').alias('id'),
                    pf.col('hp_term').alias('label')
                )).alias('diseaseModelAssociatedHumanPhenotypes')
            )
            .select('model_id', 'disease_id', 'diseaseModelAssociatedHumanPhenotypes')
        )

        self.evidence = (
            self.disease_model_summary

            # Filter out the associations with a low score. Some associations lack this score and are kept.
            .filter(~(pf.col('disease_model_max_norm') < score_cutoff))

            # Add the mouse gene mapping information. The mappings are not necessarily one to one, because a single MGI
            # can map to multiple Ensembl mouse genes. When this happens, join will handle the necessary explosions, and
            # a single row from the original table will generate multiple evidence strings.
            .join(self.mgi_gene_id_to_ensembl_mouse_gene_id, on='mgi_gene_id', how='inner')  # `targetInModel`.
            # Add the human gene mapping information. This is added in two stages: MGI → HGNC → Ensembl human gene.
            # Similarly to mouse gene mappings, at each stage there is a possibility of a row explosion.
            .join(self.mouse_gene_to_human_gene, on='mgi_gene_id', how='inner')
            .join(self.hgnc_gene_id_to_ensembl_human_gene_id, on='hgnc_id', how='inner')  # `targetFromSourceId`.
            .drop('mgi_gene_id', 'hgnc_id')

            # Add all mouse phenotypes of the model → `diseaseModelAssociatedModelPhenotypes`.
            .join(all_mouse_phenotypes, on='model_id', how='left')
            # Add the matched model/disease human phenotypes → 'diseaseModelAssociatedHumanPhenotypes`.
            .join(matched_human_phenotypes, on=['model_id', 'disease_id'], how='left')

            # Strip trailing modifiers from the model ID.
            # For example: 'MGI:6274930#hom#early' → 'MGI:6274930'.
            .withColumn(
                'biologicalModelId',
                pf.split(pf.col('model_id'), '#').getItem(0)
            )
            .drop('model_id')
            # Convert the percentage score into fraction.
            .withColumn('resourceScore', pf.col('resourceScore') / 100.0)
            # Rename the disease data columns.
            .withColumnRenamed('disease_id', 'diseaseFromSourceId')
            .withColumnRenamed('disease_term', 'diseaseFromSource')
            # Add constant value columns.
            .withColumn('datasourceId', pf.lit('phenodigm'))
            .withColumn('datatypeId', pf.lit('animal_model'))

            # Ensure stable column order.
            .select('biologicalModelAllelicComposition', 'biologicalModelGeneticBackground', 'biologicalModelId',
                    'datasourceId', 'datatypeId', 'diseaseFromSource', 'diseaseFromSourceId',
                    'diseaseModelAssociatedHumanPhenotypes', 'diseaseModelAssociatedModelPhenotypes', 'resourceScore',
                    'targetFromSourceId', 'targetInModel')
        )

    def write_evidence_strings(self, evidence_strings_filename):
        """Dump the Spark evidence dataframe into a temporary directory as separate JSON chunks. Collect and combine
        them to obtain the final output file. The order of the evidence strings is not maintained, and they are returned
        in random order as collected by Spark."""
        with tempfile.TemporaryDirectory() as tmp_dir_name, open(evidence_strings_filename, 'wb') as outfile:
            (
                self.evidence.write
                .format('json').mode('overwrite').option('compression', 'org.apache.hadoop.io.compress.GzipCodec')
                .save(tmp_dir_name)
            )
            for json_chunk_filename in [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]:
                with open(os.path.join(tmp_dir_name, json_chunk_filename), 'rb') as json_chunk:
                    shutil.copyfileobj(json_chunk, outfile)

    def process_all(self, output, score_cutoff, use_cached):
        if not use_cached:
            self.logger.info('Update the HGNC/MGI/SOLR cache.')
            self.update_cache()

        self.logger.info('Load gene mappings and SOLR data from local cache.')
        self.load_data_from_cache()

        self.logger.info('Build the evidence strings.')
        self.generate_phenodigm_evidence_strings(score_cutoff)

        self.logger.info('Collect and write the evidence strings.')
        self.write_evidence_strings(output)


def main(cache_dir, output, score_cutoff, use_cached=False, log_file=None):
    # Initialize the logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    # Process the data.
    PhenoDigm(logging, cache_dir).process_all(output, score_cutoff, use_cached)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--cache-dir', help='Directory to store the HGNC/MGI/SOLR cache files in.', required=True)
    parser.add_argument('--output', help='Name of the json.gz file to output the evidence strings into.', required=True)
    parser.add_argument('--score-cutoff', help=(
        'Discard model-disease associations with the `disease_model_max_norm` score less than this value. The score '
        'range is 0 to 100.'
    ), type=float, default=DEFAULT_ASSOCIATION_SCORE_CUTOFF)
    parser.add_argument('--use-cached', help='Use the existing cache and do not update it.', action='store_true')
    parser.add_argument('--log-file', help='Optional filename to redirect the logs into.')
    args = parser.parse_args()
    main(args.cache_dir, args.output, args.score_cutoff, args.use_cached, args.log_file)
