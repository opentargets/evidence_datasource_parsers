#!/usr/bin/env python3
"""Evidence parser for the Genomics England PanelApp data."""

import argparse
import json
import logging
import multiprocessing as mp
import re

import requests
from ontoma import OnToma
from pyspark.sql import SparkSession
from pyspark.sql.functions import (
    col, lit, when, array_distinct, split, explode, udf, regexp_extract, trim, regexp_replace, element_at, length
)
from pyspark.sql.types import StringType, ArrayType


class PanelAppEvidenceGenerator:

    REPLACE_BEFORE_SPLITTING = {
        # Fix separator errors for specific records in the source data.
        r'\(HP:0006574;\);': r'(HP:0006574);',
        r'Deafness, autosomal recessive; 12': r'Deafness, autosomal recessive, 12',
        # Fix cases like "Aarskog-Scott syndrome, 305400Mental retardation, X-linked syndromic 16, 305400", where
        # several phenotypes are glued to each other due to a formatting error.
        r'(\d{6})([A-Za-z])': r'$1;$2',
    }

    TO_REMOVE = (
        r' \(no OMIM number\)',
        r' \(NO phenotype number in OMIM\)',
        r'(No|no) OMIM (phenotype|number|entry)',
        r'[( ]*(from )?PMID:? *\d+[ ).]*'
    )

    ONTOLOGY_REGEXP = (
        r'(,? *((HP|MONDO)'
        r'[:_]?'  # Optional separator
        r'\d+),? *)'
    )

    def __init__(self):
        self.spark = SparkSession.builder.appName('panelapp_parser').getOrCreate()

    def generate_panelapp_evidence(
            self,
            input_file: str,
            output_file: str,
    ):
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

        panelapp_df.select('phenotype').coalesce(1).write.format('csv').mode('overwrite').option('compression', 'gzip').save(output_file)
        return

        # Drop records where phenotypes are empty after cleaning up.
        panelapp_df = panelapp_df.filter((length(col('phenotype')) != ''))

        panelapp_df = (
            panelapp_df

            # Extract HP/MONDO ontology codes and remove them from the phenotype string.
            .withColumn('ontology_curie', regexp_extract(col('phenotype'), r',? *((HP|MONDO)[:_]\d+),? *', 1))
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'(,? *((HP|MONDO)[:_]\d+),? *)', ''))
            # Extract OMIM code and remove it from the phenotype string.
            .withColumn('omim_code', regexp_extract(col('phenotype'), r'(\d{6})', 1))
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'OMIM(\d{6})', ''))

            # Clean phenotype from special characters.
            #
            .withColumn('phenotype', trim(regexp_replace(col('phenotype'), r',? *$', '')))

            .persist()
        )

        # Map literature mappings
        literature_mappings = PanelAppEvidenceGenerator.build_literature_mappings(
            panelapp_df.select('Panel Id').toPandas()['Panel Id'].unique()
        )
        panelapp_df = panelapp_df.withColumn(
            # 'literature', self._translate(literature_mappings, col('Panel Name'), col('Symbol'))
            'literature', PanelAppEvidenceGenerator.translate(literature_mappings)('Panel Name', 'Symbol')
        )

        # Save data
        panelapp_df.coalesce(1).write.format('json').mode('overwrite').option('compression', 'gzip').save(output_file)

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
                    gene_symbol: list({re.match(r'(\d{8})', e)[0] for e in pubs if re.match(r'(\d{8})', e) is not None})
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
    PanelAppEvidenceGenerator().generate_panelapp_evidence(
        input_file=args.input_file, output_file=args.output_file
    )



########################################################################################################################



class PanelAppEvidenceGenerator:

    def __init__(self, phenotypesMappings, limit=None):
        # Create OnToma object
        self.otmap = OnToma()

        # Create spark session
        self.spark = (
            SparkSession.builder
            .appName('evidence_builder')
            .getOrCreate()
        )
        self.dataframe = None

        # Initialize mapping variables
        self.diseaseMappings = phenotypesMappings
        self.codesMappings = None
        self.limit = limit

    def writeEvidenceFromSource(self, inputFile, skipMapping):
        '''
        Processing of the dataframe to build all the evidences from its data

        Args:
            dataframe (pyspark.DataFrame): Initial .tsv file
        Returns:
            evidences (array): Object with all the generated evidences strings
        '''

        # Reading and filtering input file
        self.dataframe = (
            self.spark.read.csv(inputFile, sep=r'\t', header=True)
            .filter(
                ((col('List') == 'green') | (col('List') == 'amber'))
                & (col('Panel Version') > 1) & (col('Panel Status') == 'PUBLIC')
            )
        )

        # Applying limit is present:
        if self.limit is not None:
            self.dataframe = self.dataframe.sample(False, 1.0, 829348).limit(self.limit)

        logging.info('Fetching publications from the API...')
        pdf = PanelAppEvidenceGenerator.buildPublications(self.dataframe.toPandas())  # TODO: write in pyspark
        pdf.dropna(axis=1, how='all', inplace=True)
        self.dataframe = self.spark.createDataFrame(pdf)
        logging.info('Publications loaded.')

        # Cleaning the phenotype related data of the dataframe
        self.dataframe = PanelAppEvidenceGenerator.cleanDataframe(self.dataframe)

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

    @staticmethod
    def buildPublications(pdf):
        '''
        Populates a dataframe with the publications fetched from the PanelApp API and cleans them to match PubMed IDs.

        Args:
            dataframe (pandas.DataFrame): DataFrame with transformed PanelApp data
        Returns:
            dataframe (pandas.DataFrame): DataFrame with an 'publications' column added
        '''
        populated_groups = []

        for (PanelId), group in pdf.groupby('Panel Id'):
            request = PanelAppEvidenceGenerator.publicationsFromPanel(PanelId)
            group['publications'] = group.apply(
                lambda X: PanelAppEvidenceGenerator.publicationFromSymbol(X.Symbol, request), axis=1
            )
            populated_groups.append(group)

        pdf = pd.concat(populated_groups, ignore_index=True, sort=False)

        cleaned_publication = []
        for row in pdf['publications'].to_list():
            try:
                tmp = [re.match(r'(\d{8})', e)[0] for e in row]
                cleaned_publication.append(list(set(tmp)))  # Removing duplicated publications
            except Exception as e:
                cleaned_publication.append([])
                continue

        pdf['publications'] = cleaned_publication

        return pdf

    @staticmethod
    def publicationsFromPanel(panelId):
        '''
        Queries the PanelApp API to obtain a list of the publications for every gene within a panelId

        Args:
            panelId (str): Panel ID extracted from the 'Panel Id' column
        Returns:
            response (dict): Response of the API containing all genes related to a panel and their publications
        '''
        try:
            url = f'http://panelapp.genomicsengland.co.uk/api/v1/panels/{panelId}/'
            res = requests.get(url).json()
            return res['genes']
        except Exception as e:
            logging.error('Query of the PanelApp API has failed.')
            return None

    @staticmethod
    def publicationFromSymbol(symbol, response):
        '''
        Returns the list of publications for a given symbol in a PanelApp query response.
        Args:
            symbol (str): Gene symbol extracted from the 'Symbol' column
            response (dict): Response of the API containing all genes related to a panel and their publications
        Returns:
            publication (list): Array with all publications for a particular gene in the corresponding Panel ID
        '''
        for gene in response:
            if gene['gene_data']['gene_symbol'] == symbol:
                publication = gene['publications']
                return publication

    @staticmethod
    def cleanDataframe(dataframe):
        '''
        Args:
            dataframe (pyspark.DataFrame): Initial .tsv data
        Returns:
            dataframe (pyspark.DataFrame): Transformed initial dataframe
        '''

        dataframe = (
            dataframe
            .withColumn(
                # NaNs and 'No OMIM phenotype' in 'Phenotypes' column --> Assignment of Panel Name
                'Phenotypes',
                when(col('Phenotypes') == 'No OMIM phenotype', col('Panel Name'))
                .when(col('Phenotypes').isNull(), col('Panel Name'))
                .otherwise(col('Phenotypes'))
            )
            # cohortPhenotypes --> array of the original string separated by phenotypes
            .withColumn('cohortPhenotypes', array_distinct(split(col('Phenotypes'), ';')))
            # phenotype --> explosion of cohortPhenotypes
            .withColumn('phenotype', explode(col('cohortPhenotypes')))
            # omimCode --> OMIM code present in 'phenotype'
            .withColumn('omimCode', regexp_extract(col('phenotype'), r'(\d{6})', 1))
            # removal of the OMIM code in 'phenotype'
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'(\d{6})', ''))
            # deleting special characters in 'phenotype'
            .withColumn('phenotype', regexp_replace(col('phenotype'), r'[^0-9a-zA-Z *]', ''))
            .withColumn('phenotype', trim(col('phenotype')))
        )

        return dataframe

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
