#!/usr/bin/env python3
"""Evidence parser for the Genomics England PanelApp data."""

import argparse
import logging
import os
import re
import tempfile

import requests
from pyspark.sql import SparkSession
from pyspark.sql.functions import (
    col, collect_set, concat, lit, when, array_distinct, split, explode, regexp_extract, trim, regexp_replace
)


class PanelAppEvidenceGenerator:

    REPLACE_BEFORE_SPLITTING = {
        # Fix separator errors for specific records in the source data.
        r'\(HP:0006574;\);': r'(HP:0006574);',
        r'Deafness, autosomal recessive; 12': r'Deafness, autosomal recessive, 12',
        r'Waardenburg syndrome, type; 3': r'Waardenburg syndrome, type 3',
        r'Ectrodactyly, ectodermal dysplasia, and cleft lip/palate syndrome; 3':
            r'Ectrodactyly, ectodermal dysplasia, and cleft lip/palate syndrome, 3',

        # Remove curly braces. They are *sometimes* (not consistently) used to separate disease name and OMIM code, for
        # example: "{Bladder cancer, somatic}, 109800".
        r'[{}]': r'',

        # Fix cases like "Aarskog-Scott syndrome, 305400Mental retardation, X-linked syndromic 16, 305400", where
        # several phenotypes are glued to each other due to formatting errors.
        r'(\d{6})([A-Za-z])': r'$1;$2',

        # Replace all tabs with spaces and multiple spaces with one space.
        r'\t': r' ',
        r' +': r' ',
    }

    TO_REMOVE = (
        r' \(no OMIM number\)',
        r' \(NO phenotype number in OMIM\)',
        r'(No|no) OMIM (phenotype|number|entry)',
        r'[( ]*(from )?PMID:? *\d+[ ).]*'
    )

    # Regular expressions for extracting ontology information.
    ONT_LEADING = r'[ ,-]*'
    ONT_SEPARATOR = r'[:_ #]*'
    ONT_TRAILING = r'[:.]*'
    OMIM_REGEXP = ONT_LEADING + r'(OMIM|MIM)?' + ONT_SEPARATOR + r'(\d{6})' + ONT_TRAILING
    OTHER_REGEXP = ONT_LEADING + r'(OrphaNet: ORPHA|Orphanet|ORPHA|HP|MONDO)' + ONT_SEPARATOR + r'(\d+)' + ONT_TRAILING

    # Regular expressions for extracting publication information.
    PMID_REGEXPS = [
        (
            r'^'        # Start of the string
            r'[\d, ]+'  # Followed by a sequence of digits, commas and spaces
            r'(?: |$)'  # Ending either with a space, or with the end of the string
        ),
        (
            r'(?:PubMed|PMID)'  # PubMed or a PMID prefix
            r'[: ]*'            # An optional separator (spaces/colons)
            r'[\d, ]'           # A sequence of digits, commas and spaces
        )
    ]

    def __init__(self):
        self.spark = SparkSession.builder.appName('panelapp_parser').getOrCreate()

    def generate_panelapp_evidence(
        self,
        input_file: str,
        output_file: str,
    ) -> None:
        # Filter and extract the necessary columns.
        panelapp_df = (
            self.spark.read.csv(input_file, sep=r'\t', header=True)
            .filter(
                ((col('List') == 'green') | (col('List') == 'amber')) &
                (col('Panel Version') > 1) &
                (col('Panel Status') == 'PUBLIC')
            )
            .select('Symbol', 'Panel Id', 'Panel Name', 'Mode of inheritance', 'Phenotypes')
        )

        # Fix typos and formatting errors which interfere with phenotype splitting.
        for regexp, replacement in self.REPLACE_BEFORE_SPLITTING.items():
            panelapp_df = panelapp_df.withColumn('Phenotypes', regexp_replace(col('Phenotypes'), regexp, replacement))

        # Split and explode the phenotypes.
        panelapp_df = (
            panelapp_df
            .withColumn('cohortPhenotypes', array_distinct(split(col('Phenotypes'), ';')))
            .withColumn('phenotype', explode(col('cohortPhenotypes')))
            .drop('Phenotypes')
        )

        # Remove specific patterns and phrases which will interfere with ontology extraction and mapping.
        for regexp in self.TO_REMOVE:
            panelapp_df = panelapp_df.withColumn('phenotype', regexp_replace(col('phenotype'), f'({regexp})', ''))

        # Extract ontology information, clean up and filter the split phenotypes.
        panelapp_df = (
            panelapp_df

            # Extract Orphanet/MONDO/HP ontology identifiers and remove them from the phenotype string.
            .withColumn('ontology_namespace', regexp_extract(col('phenotype'), self.OTHER_REGEXP, 1))
            .withColumn('ontology_namespace', regexp_replace(col('ontology_namespace'), 'OrphaNet: ORPHA', 'ORPHA'))
            .withColumn('ontology_id', regexp_extract(col('phenotype'), self.OTHER_REGEXP, 2))
            .withColumn(
                'ontology',
                when(
                    col('ontology_namespace').isNotNull() & col('ontology_id').isNotNull(),
                    concat(col('ontology_namespace'), lit(':'), col('ontology_id'))
                )
            )
            .withColumn('phenotype', regexp_replace(col('phenotype'), f'({self.OTHER_REGEXP})', ''))

            # Extract OMIM identifiers and remove them from the phenotype string.
            .withColumn('omim_id', regexp_extract(col('phenotype'), self.OMIM_REGEXP, 2))
            .withColumn('omim', when(col('omim_id').isNotNull(), concat(lit('OMIM:'), col('omim_id'))))
            .withColumn('phenotype', regexp_replace(col('phenotype'), f'({self.OMIM_REGEXP})', ''))

            # Choose one of the ontology identifiers, keeping OMIM as a priority.
            .withColumn('diseaseFromSourceId', when(col('omim').isNotNull(), col('omim')).otherwise(col('ontology')))
            .drop('ontology_namespace', 'ontology_id', 'ontology', 'omim_id', 'omim')

            # Clean up the final split phenotypes.
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'\(\)', ''))
            .withColumn('phenotype', trim(col('phenotype')))
            .withColumn('phenotype', when(col('phenotype') != '', col('phenotype')))

            # Remove empty or low quality records.
            .filter(
                ~(
                    # There is neither a phenotype string nor an ontology identifier.
                    ((col('phenotype').isNull()) & (col('diseaseFromSourceId').isNull())) |
                    # The name of the phenotype string starts with a question mark, indicating low quality.
                    col('phenotype').startswith('?')
                )
            )
            .persist()
        )

        # Fetch and join literature references.
        all_panel_ids = panelapp_df.select('Panel Id').toPandas()['Panel Id'].unique()
        literature_references = self.fetch_literature_references(all_panel_ids)
        panelapp_df = panelapp_df.join(literature_references, on=['Panel Id', 'Symbol'], how='left')

        # Populate final evidence string structure.
        panelapp_df = (
            panelapp_df
            .withColumnRenamed('Mode of inheritance', 'allelicRequirements')
            .withColumnRenamed('List', 'confidence')
            .withColumn('datasourceId', lit('genomics_england'))
            .withColumn('datatypeId', lit('genetic_literature'))
            .withColumnRenamed('phenotype', 'diseaseFromSource')
            # diseaseFromSourceId populated above
            # literature populated above
            .withColumnRenamed('Panel Id', 'studyId')
            .withColumnRenamed('Panel Name', 'studyOverview')
            .withColumnRenamed('Symbol', 'targetFromSourceId')
        )

        # Save data.
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            (
                panelapp_df.coalesce(1).write.format('json').mode('overwrite')
                .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
            )
            json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
            assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
            os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)

    def fetch_literature_references(self, all_panel_ids):
        """Queries the PanelApp API to extract all literature references for (panel ID, gene symbol) combinations."""
        publications = []  # Contains tuples of (panel ID, gene symbol, PubMed ID).
        for panel_id in all_panel_ids:
            logging.info(f'Fetching literature references for panel {panel_id!r}.')
            url = f'https://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}'
            for gene in requests.get(url).json()['genes']:
                for publication_string in gene['publications']:
                    publications.extend([
                        (panel_id, gene['gene_data']['gene_symbol'], pubmed_id)
                        for pubmed_id in self.extract_pubmed_ids(publication_string)
                    ])
        # Group by (panel ID, gene symbol) pairs and convert into a PySpark dataframe
        return (
            self.spark
            .createDataFrame(publications, schema=('Panel ID', 'Symbol', 'literature'))
            .groupby(['Panel ID', 'Symbol'])
            .agg(collect_set('literature').alias('literature'))
        )

    def extract_pubmed_ids(self, publication_string):
        """Parses the publication information from the PanelApp API and extracts PubMed IDs."""
        publication_string = publication_string.strip()
        pubmed_ids = []

        for regexp in self.PMID_REGEXPS:  # For every known representation pattern...
            for occurrence in re.findall(regexp, publication_string):  # For every occurrence of this pattern, if any...
                pubmed_ids.extend(re.findall(r'(\d+)', occurrence))  # Extract all digit sequences = PMIDs.

        # Filter out '0' as a value, because it is a placeholder for a missing ID.
        return {pubmed_id for pubmed_id in pubmed_ids if pubmed_id != '0'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--input-file', required=True, type=str,
        help='Input TSV file with the table containing association details.'
    )
    parser.add_argument(
        '--output-file', required=True, type=str,
        help='Output JSON file (gzip-compressed) with the evidence strings.'
    )
    args = parser.parse_args()
    PanelAppEvidenceGenerator().generate_panelapp_evidence(input_file=args.input_file, output_file=args.output_file)
