#!/usr/bin/env python3

import argparse
import logging
import gzip
import json
import multiprocessing as mp
import re
from sys import stderr

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.functions import (
    col, lit, when, array_distinct, split, explode, udf, regexp_extract, trim, regexp_replace, element_at
)
from pyspark.sql.types import StringType, ArrayType
import requests

from ontoma import OnToma

class PanelAppEvidenceGenerator():
    def __init__(self):
        # Initialize spark session
        local=True
        if local:
            sparkConf = (
                SparkConf()
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
        else:
            sparkConf = (
                SparkConf()
                .set('spark.driver.maxResultSize', '0')
                .set('spark.debug.maxToStringFields', '2000')
                .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
            )
            self.spark = (
                SparkSession.builder
                .config(conf=sparkConf)
                .getOrCreate()
        )
        self.evidence = None

    def generate_panelapp_evidence(
        self,
        input_file: str,
        output_file: str,
    ):
        panelapp_df = (
            self.spark.read.csv(input_file, sep=r'\t', header=True)
            .filter(
                ((col('List') == 'green') | (col('List') == 'amber')) &
                (col('Panel Version') > 1) &
                (col('Panel Status') == 'PUBLIC')
            )
            .select(
                'Symbol', 'Panel Id', 'Panel Name',
                'Mode of inheritance', 'Phenotypes'
            )
            # If Phenotypes is empty, Panel Name is assigned
            .withColumn(
                'Phenotypes',
                when(((col('Phenotypes') == 'No OMIM phenotype') | (col('Phenotypes').isNull())), col('Panel Name'))
                .otherwise(col('Phenotypes'))
            )
            .withColumn('cohortPhenotypes', array_distinct(split(col('Phenotypes'), ';')))
            .withColumn('phenotype', explode(col('cohortPhenotypes')))
            # Extract OMIM code
            .withColumn('omim_code', regexp_extract(col('phenotype'), r'(\d{6})', 1))
            # Remove OMIM code from phenotype
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'(\d{6})', ''))
            # Clean phenotype from special characters
            .withColumn('phenotype',
                trim(
                    regexp_replace(col('phenotype'), r'[^0-9a-zA-Z -]', '')
            ))
            .drop('Phenotypes')
            .persist()
        )

        # Map literature mappings
        literature_mappings = PanelAppEvidenceGenerator.build_literature_mappings(
            panelapp_df.select('Panel Id').toPandas()['Panel Id'].unique()
        )
        print(list(literature_mappings.items())[:4])
        panelapp_df = panelapp_df.withColumn(
            #'literature', self._translate(literature_mappings, col('Panel Name'), col('Symbol'))
            'literature', PanelAppEvidenceGenerator.translate(literature_mappings)('Panel Name', 'Symbol')
        )

        # Save data
        panelapp_df.write.format('json').mode('overwrite').option('compression', 'gzip').save(output_file)

    @staticmethod
    def build_literature_mappings(
        panel_ids: 'list[str]'
    ) -> dict:
        literature_mappings = dict()
        for panel_id in panel_ids:
            res = PanelAppEvidenceGenerator.query_literature(panel_id)
            for gene in res:
                gene_symbol = gene['gene_data']['gene_symbol']
                pubs = gene.get('publications', default=None)
                literature_mappings[panel_id] = {
                    gene_symbol : list({re.match(r'(\d{8})', e)[0] for e in pubs if re.match(r'(\d{8})', e) is not None})
                }
        return literature_mappings

    @staticmethod
    def query_literature(
        panel_id: str
    ) -> 'list[str]':
        """
        Queries the PanelApp API to obtain a list of the publications for every gene within a panelId
        """
        try:
            url = f'http://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/'
            return requests.get(url).json()['genes']
        except Exception:
            print('Query of the PanelApp API has failed.')
    
    @staticmethod
    def translate(
        literature_mappings: dict
    ):
        """
        Mapping panel/gene pairs to literature
        """
        def translate_(panel, gene):
            return literature_mappings.get(panel).get(gene)
        return udf(translate_, ArrayType(StringType()))
    
    @staticmethod
    @udf(ArrayType(StringType()))
    def _translate(self, literature_mappings, panel, gene):
            return literature_mappings.get(panel).get(gene)


PanelAppEvidenceGenerator().generate_panelapp_evidence(input_file="Some_genes.tsv", output_file='output/first_iter.json.gz')