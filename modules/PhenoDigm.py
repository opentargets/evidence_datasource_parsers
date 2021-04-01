import os
import gzip
import json
import argparse
import logging
import pandas as pd

import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from pyspark import SparkContext
import re




class IPMC_solr_parser(object):
    '''
    This class retrieves data from the IPMC solr API.
    
    * Returns array of documents.
    * If target folder is specified, json files are saved into that folder.
    '''
    
    ipmc_solr_host = 'http://www.ebi.ac.uk/mi/impc/solr/phenodigm/select'

    
    def __init__(self, target_folder, rows=20000, limit=None):
        """
        Storing basic values when initializing object
        
        Args:
        rows (int): Number of solr documents returned in a single query
        limit (int): Maximum number of returned document. If None, all documents are returned
        target_folder (string): Folder into which the data is saved.
        
        Returns: 
        None
        """
        
        self.rows = rows
        self.limit = limit
        self.target_folder = target_folder
        
        
    def fetch_data(self, data_type=None):
        """
        Fetching data to the specified location. If data type is not specified, 
        all types are retrieved and the files are saved directly to the root folder
        
        Args:
        data_type (string): data type to return match .type == 'data_type'
        
        Returns:
        None
        """
        
        # Based on the data type, we update the output folder:
        data_folder = f'{self.target_folder}/type.{data_type}' if data_type else self.target_folder
            
        # Create folder:
        os.makedirs(data_folder, exist_ok=True) 
        
        # Initialize counter:
        start = 0
        total = 0
        numFound = 1
        chunk = 0
        limit = None

        logging.info(f'Retrieving data from IPMC solr: {self.ipmc_solr_host}')
        logging.info(f'Retrieving {self.rows} documents at a time.')
        logging.info(f'Specified data type: {data_type}')
            
        # Retrieving data from IMPC:
        while True:
            
            # Retrieve data:
            solr_data = self.query_solr(start, data_type)
       
            # If limit is not set, return all data:
            if not limit:
                limit = solr_data['response']['numFound'] if not self.limit else self.limit
            
            # If we don't find any items, we break:
            if solr_data['response']['numFound'] == 0:
                break
                
            # Write data to file:
            with gzip.open(f'{data_folder}/IMPC_solr_dump.{chunk:03}.json.gz', 'tw') as f:
                for doc in solr_data['response']['docs']:
                    json.dump(doc, f)
                    f.write('\n')

            # Incrementing counters:
            start += self.rows
            chunk += 1
            total += len(solr_data['response']['docs'])
            
            # If the length of the dataframe reaches the limit, we exit:
            if (limit <= total) or (total == solr_data['response']['numFound']):
                break
            
            # Log progress
            logging.debug(f'Chunk {chunk} done. Number of retrieved documents: {total}.')

        logging.info(f'Retrieval finished. Number of documents: {total}')
            

    # Use @retry decorator to ensure that errors like the query failing because server was overloaded, are handled correctly and the request is retried
    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def query_solr(self, start, data_type=None):

        # Building request:
        params = {'q': '*:*','start': start,'rows': self.rows}

        if data_type:
            params['fq'] = f'type:{data_type}'
            
        # Query
        response = requests.get(self.ipmc_solr_host, params=params, timeout=30)

        # Check for erroneous HTTP response statuses
        response.raise_for_status()
        response_data = response.json()
        return response_data


def get_solr_data(target_folder):
    """
    This function retrieves solr data to a defined location.

    Parameters:
    target_folder (str) - location into witch solr files are saved

    Returns:
    None
    """

    # List of data types to retrieve:
    data_types = ['gene','gene_gene','mouse_model','disease_model_summary',
                  'disease_gene_summary','disease','ontology_ontology','ontology']

    # Initialize impc solr object:
    impc_solr_retriever = IPMC_solr_parser(target_folder=target_folder)

    # Looping through all 
    for data_type in data_types:
        logging.info(f'fetching: {data_type}')
        impc_solr_retriever.fetch_data(data_type=data_type)


def process_solr_data(solr_data_folder, output_file):
    """
    This function reads the downloaded solr documents into spark dataframes
    Does the required actions - filtering, joining, aggregation
    And saves the output file.

    Parameters:
    solr_data_folder (str): frolder from which spark reads data
    output_file (str): output file name

    Returns:
    None
    """

    global spark

    # Creating spark session:
    spark = (pyspark.sql.SparkSession
        .builder
        .appName("phenodigm_parser")
        .getOrCreate()
    )


    logging.info('Spark version: ', spark.version)


    """
    The steps below 
     * read the gene set
     * filter for human genes (rows with hgnc_id)
     * join with gene-gene linking table

    Resulting table has 
     * human HGNC gene IDs
     * human HGNC gene symbol
     * mouse MGI gene id
    """
    logging.info('Reading gene table...')

    human_genes = (
        spark.read.json(f'{solr_data_folder}/type.gene')
        .select('hgnc_gene_id', 'hgnc_gene_symbol')
        .filter(col('hgnc_gene_id').isNotNull())
    )

    genes_table = (
        spark.read.json(f'{solr_data_folder}/type.gene_gene')
        .select('hgnc_gene_id', 'gene_id')
        .join(human_genes, on='hgnc_gene_id', how='inner')
    )

    logging.info(f'Number of human genes with mouse orthologue: { genes_table.select('hgnc_id').distinct().count()}')

    """
    Ontology table contains the mapping between human an mouse
    phenotypes. Only those mouse phenotypes included that have human correspondent

    The mouse phenotype term is not included - that value comes from the models table
    """
    ontolgy_table = (
        spark.read.json(f'{solr_data_folder}/type.ontology_ontology')
        .select('hp_id','hp_term','mp_id')
    )

    logging.info(f"Number of human phenotypes: {ontolgy_table.select('hp_id').distinct().count()}")
    logging.info(f"Number of mouse phenotypes: {ontolgy_table.select('mp_id').distinct().count()}")
    logging.info(f"Number of human to mouse phenotype mappings: {ontolgy_table.count()}")

    ##
    ## Processing mouse models:
    ##

    mouse_model_table = (
        # Read mouse model:
        spark.read.json(f'{solr_data_folder}/type.mouse_model')

        # Select columns, rename:
        .select('model_id','model_phenotypes', 'marker_id')
        .withColumnRenamed('marker_id', 'gene_id')

        # The phenotypes column exploded and the id/term pair extracted:
        .withColumn('model_phenotype', explode(col('model_phenotypes')))
        .withColumn('parsed_phenotype', parse_phenotypes(col('model_phenotype')))
        .drop('model_phenotypes')
        .select('model_id', 'gene_id', 'parsed_phenotype.mp_id', 'parsed_phenotype.mp_term')

        # Joining ontology table that maps mouse phenotypes to human phenotypes:
        .join(ontolgy_table, on='mp_id', how='left')

        # Dropping phenotypes without human mappings:
        .filter(col('hp_id').isNotNull())
    )

    # Aggregating phenotypes and joining genes:
    mouse_model_phenotype_parsed = (
        mouse_model_table

        # Grouping rows by model id and mouse gene id (human genes are not there yet)
        .groupby('model_id','gene_id')
        .agg(

            # Mouse phenotypes are pooled into a list of struct:
            collect_set(struct(           
                col("mp_id").alias('id'), 
                col('mp_term').alias('label')
            )).alias('diseaseModelAssociatedModelPhenotypes'),

            # Human phenotypes are pooled into a list of struct:
            collect_set(struct(           
                col("hp_id").alias('id'), 
                col('hp_term').alias('label')
            )).alias('diseaseModelAssociatedHumanPhenotypes')        
        )

        # Joining with human genes:
        .join(genes_table, on='gene_id', how='inner')
    )

    ##
    ## Opening disease table and join with models:
    ##
    disease_model_table = (
        spark.read.json(f'{solr_data_folder}/type.disease_model_summary/')
        .drop(*['association_curated', 'marker_locus', 'marker_symbol','model_source','type', 'disease_model_avg_norm', 'marker_id', 'marker_num_models'])
        .join(mouse_model_phenotype_parsed, on='model_id', how='inner')
    )


def main():

    # Parsing arguments:
    parser = argparse.ArgumentParser(description='Evidence parser for animal models sources from PhenoDigm')

    # Two of the flags are mutually exclusive:
    parser.add_argument("-u", '--use-cache', help='To use cached data instead of fetching from IUPMC solr', 
        action="store_true", dest="update_cache", default=False)
    parser.add_argument("-c", '--cache-dir', help='Location where the cached datafiles are saved/read.',
        type=str, dest="cache_dir", required=True)
    parser.add_argument("-o", '--output_file', help='Name of the output file (gzipped json).',
        type=str, dest="output_file", required=True)

    # log file is optional:
    parser.add_argument("-l", '--logFile', help='Optional filename for logfile.', required=False)
    args = parser.parse_args()

    # Initialize logger:
    # Initialize logger based on the provided logfile. 
    # If no logfile is specified, logs are written to stderr 
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    
    cache_dir = args.cache_dir
    output_file = args.output_file

    # Retrieving samples from every data types:
    if args.update_cache:
        get_solr_data(cache_dir)

    # Processing data and generate evidence:
    process_solr_data(cache_dir, output_file)
    


if __name__ == '__main__':
    main()

