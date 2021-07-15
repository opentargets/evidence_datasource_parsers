#!/usr/bin/env python3
"""Evidence parser for the Genomics England PanelApp data."""

import argparse
import json
import logging
import multiprocessing as mp
import os
import re
import tempfile

import requests
from ontoma import OnToma
from pyspark.sql import SparkSession
from pyspark.sql.functions import (
    col, collect_set, concat, lit, when, array_distinct, split, explode, udf, regexp_extract, trim, regexp_replace, element_at
)
from pyspark.sql.types import StringType, ArrayType


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
            .withColumn(
                'ontology',
                concat(
                    regexp_extract(col('phenotype'), self.OTHER_REGEXP, 1),
                    lit(':'),
                    regexp_extract(col('phenotype'), self.OTHER_REGEXP, 2)
                )
            )
            .withColumn('ontology', regexp_replace(col('ontology'), 'OrphaNet: ORPHA', 'ORPHA'))
            .withColumn('ontology', regexp_replace(col('ontology'), '^:$', ''))
            .withColumn('phenotype', regexp_replace(col('phenotype'), f'({self.OTHER_REGEXP})', ''))

            # Extract OMIM identifiers and remove them from the phenotype string.
            .withColumn(
                'omim_id',
                concat(
                    lit('OMIM:'),
                    regexp_extract(col('phenotype'), self.OMIM_REGEXP, 2)
                )
            )
            .withColumn('omim_id', regexp_replace(col('omim_id'), '^OMIM:$', ''))
            .withColumn('phenotype', regexp_replace(col('phenotype'), f'({self.OMIM_REGEXP})', ''))

            # Choose one of the ontology identifiers, keeping OMIM as a priority.
            .withColumn('diseaseFromSourceId', when(col('omim_id') != '', col('omim_id')).otherwise(col('ontology')))
            .drop('ontology', 'omim_id')

            # Clean up the final split phenotypes.
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'\(\)', ''))
            .withColumn('phenotype', trim(col('phenotype')))

            # Remove empty or low quality records.
            .filter(
                ~(
                    # There is neither a string, nor an OMIM identifier, nor some other ontology identifier.
                    ((col('phenotype') == '') & (col('omim_id') == '') & (col('ontology') == '')) |
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



########################################################################################################################



class PanelAppEvidenceGenerator:

    def writeEvidenceFromSource(self, inputFile, skipMapping):
        # Map the diseases to an EFO term if necessary
        if not skipMapping:
            self.dataframe = self.diseaseMappingStep()
        else:
            logging.info('Disease mapping has been skipped.')
            self.dataframe = self.dataframe.withColumn(
                'ontomaUrl',
                lit(None)
            )

        logging.info('Generating evidence:')
        evidences = (
            self.dataframe
            # Removing redundant evidence after the explosion of phenotypes
            .dropDuplicates(['Panel Id', 'Symbol', 'ontomaUrl', 'cohortPhenotypes'])
            # Build evidence strings per row
            .rdd.map(PanelAppEvidenceGenerator.parseEvidenceString)
            .collect()  # list of dictionaries
        )

        for evidence in evidences:
            # Delete empty keys
            if skipMapping or evidence['diseaseFromSourceMappedId'] is None:
                del evidence['diseaseFromSourceMappedId']

            if not evidence['diseaseFromSourceId']:
                del evidence['diseaseFromSourceId']

            if len(evidence['literature']) == 0:
                del evidence['literature']

            if evidence['allelicRequirements'] == [None]:
                del evidence['allelicRequirements']

        return evidences

    def diseaseMappingStep(self):
        '''
        Runs the disease mapping step from phenotypes and OMIM codes to an EFO ID - 3 steps:
            1. Querying OnToma with all distinct codes and phenotypes
            2. Xref between the results of every phenotype/OMIM code pair for a better coverage
            3. Mappings are built into the dataframe

        Args:
            dataframe (pyspark.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pyspark.DataFrame): DataFrame with the mapping results filtered by only matches
        '''

        omimCodesDistinct = list(self.dataframe.select('omimCode').distinct().toPandas()['omimCode'])
        phenotypesDistinct = list(self.dataframe.select('phenotype').distinct().toPandas()['phenotype'])

        if self.diseaseMappings is None:
            # Checks whether the dictionary is not provided as a parameter
            try:
                pool = mp.Pool(mp.cpu_count())
                self.diseaseMappings = pool.map(self.diseaseToEfo, phenotypesDistinct)  # list of dictionaries
                self.diseaseMappings = {k: v for dct in self.diseaseMappings for (k, v) in self.diseaseMappings.items()}
                pool.close()
            except Exception as e:
                logging.error(f'Processing the mappings without multithreading. Following error occurred: \n{e}.')
                self.diseaseMappings = self.diseaseToEfo(*phenotypesDistinct)
        else:
            logging.info('Disease mappings have been imported.')

        self.codesMappings = self.diseaseToEfo(*omimCodesDistinct, dictExport='codesToEfo_results.json')
        (
            self.dataframe
            .select('omimCode', 'phenotype')
            .toPandas()
            .drop_duplicates()
            .apply(
                lambda row: self.phenotypeCodePairCheck(row['omimCode'], row['phenotype']), axis=1
            )
        )
        logging.info('Disease mappings have been checked.')

        # Add new columns: ontomaResult, ontomaUrl, ontomaLabel
        phenotypesMappings = self.diseaseMappings
        udfBuildMapping = udf(
            lambda X: PanelAppEvidenceGenerator.buildMapping(X, phenotypesMappings),
            ArrayType(StringType()))
        self.dataframe = (
            self.dataframe
            .withColumn('ontomaResult', udfBuildMapping(col('phenotype'))[0])
            .withColumn('ontomaUrl', element_at(split(udfBuildMapping(col('phenotype'))[1], '/'), -1))
            .withColumn('ontomaLabel', udfBuildMapping(col('phenotype'))[1])
        )

        return self.dataframe

    def diseaseToEfo(self, *iterable, dictExport='diseaseToEfo_results.json'):
        '''
        Queries the OnToma utility to map a phenotype to an EFO ID.

        Args:
            iterable (array): Array or column of a dataframe containing the strings to query
            dictExport (str): Name of the output file where the OnToma queries will be saved
        Returns:
            mappings (dict): Output file. Keys: queried term (phenotype or OMIM code), Values: OnToma output
        '''
        mappings = dict()
        for e in iterable:
            try:
                ontomaResult = self.otmap.find_term(e, verbose=True)
                if ontomaResult is not None:
                    mappings[e] = ontomaResult
                else:
                    mappings[e] = {
                        'term': None,
                        'label': None,
                        'source': None,
                        'quality': None,
                        'action': None
                    }
            except Exception as error:
                logging.error(f'{e} mapping has failed.')
                continue

        with open(dictExport, 'w') as outfile:
            json.dump(mappings, outfile)

        return mappings

    def phenotypeCodePairCheck(self, phenotype, omimCode):
        '''
        Among the Fuzzy results of a phenotype query, it checks if the phenotype and the respective code points
        to the same EFO term and changes the dictionary with the mappings

        Args:
            phenotype (str): Phenotype of each phenotype-OMIM code pair
            omimCode (str): Code of each phenotype-OMIM code pair
        '''
        try:
            if phenotype == '':
                return
            phenotypeResult = self.diseaseMappings[phenotype]
            if phenotypeResult is None:
                return

            if phenotypeResult['quality'] == 'fuzzy':
                codeResult = self.codesMappings[omimCode]
                if codeResult is None:
                    return

                if codeResult['term'] == phenotypeResult['term']:
                    # If both EFO terms coincide, the phenotype in mappingsDict becomes a match
                    self.diseaseMappings[phenotype]['quality'] = 'match'
                    self.diseaseMappings[phenotype]['action'] = 'checked'
        except Exception as e:
            logging.error(f'No OMIM code for phenotype: {phenotype}')

    def buildMapping(phenotype, phenotypesMappings):
        if phenotype in phenotypesMappings.keys():
            ontomaResult = phenotypesMappings[phenotype]['quality']
            ontomaUrl = phenotypesMappings[phenotype]['term']
            ontomaLabel = phenotypesMappings[phenotype]['label']
            # TO-DO: https://stackoverflow.com/questions/42980704/pyspark-create-new-column-with-mapping-from-a-dict
            return ontomaResult, ontomaUrl, ontomaLabel

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                'datasourceId': 'genomics_england',
                'datatypeId': 'genetic_literature',
                'confidence': row['List'],
                'diseaseFromSource': row['phenotype'],
                'diseaseFromSourceId': row['omimCode'],
                'diseaseFromSourceMappedId': row['ontomaUrl'],
                'cohortPhenotypes': row['cohortPhenotypes'],
                'targetFromSourceId': row['Symbol'],
                'allelicRequirements': [row['Mode of inheritance']],
                'studyId': row['Panel Id'],
                'studyOverview': row['Panel Name'],
                'literature': row['publications']
            }
            return evidence
        except Exception as e:
            logging.error(f'Evidence generation failed for row: {row}')
            raise
