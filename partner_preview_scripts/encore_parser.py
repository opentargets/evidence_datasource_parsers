#!/usr/bin/env python3
"""Parser for data submitted from the Encore team."""

import argparse
import json
import logging
import os
import sys
from functools import reduce


from pyspark.sql import SparkSession
from pyspark.sql.functions import col, struct, regexp_replace, udf, lit, expr, split, explode, array
from pyspark.sql.types import StringType, StructType, StructField, ArrayType, FloatType
from pyspark.sql.dataframe import DataFrame

from common.evidence import initialize_sparksession, write_evidence_strings, GenerateDiseaseCellLines

class EncoreEvidenceGenerator:
    """This parser generates disease/target evidence based on datafiles submitted by the Encore team.

    Input for the parser:
        - JSON configuration describing the experiments.

    For each experiment the following input files are expected:
        - lfc_file: log fold change values for every gene pairs, for every experiment.
        - bliss_file: the analysed cooperative effect for each gene paris calculated by the BLISS analysis.
        - gemini_file: the analysed cooperative effect for each gene paris calculated by the Gemini analysis.

    The parser joins the lfc data with the cooperativity measurements by gene pairs and cell line.
    Apply filter on gene pairs:
        - only extract gene pairs with significant log fold change.
        - only extract gene pairs with significant cooperative effect.
    """
    def __init__(self, spark: SparkSession, cell_passport_df: DataFrame, shared_parameters: dict) -> None:
        self.spark = spark
        self.cell_passport_df = cell_passport_df

        # Parsing paramters:
        self.log_fold_change_cutoff_p_val = shared_parameters["logFoldChangeCutoffPVal"]
        self.log_fold_change_cutoff_FDR = shared_parameters["logFoldChangeCutoffFDR"]
        self.interaction_cutoff_p_val = shared_parameters["interactionCutoffPVal"]
        self.data_folder = shared_parameters["data_folder"]
        self.release_version = shared_parameters["releaseVersion"]
        self.release_date = shared_parameters["releaseDate"]

    @staticmethod
    @udf(ArrayType(StructType([
        StructField("targetFromSourceId", StringType(), False),
        StructField("targetRole", StringType(), False),
        StructField("interactingTargetFromSourceId", StringType(), False),
        StructField("interactingTargetRole", StringType(), False),
    ])))
    def parse_targets(gene_pair: str, gene_role: str) -> dict:
        """The gene pair string is split and assigned to the relevant role + exploding into evidence of two targets.

        gene pair: 'SHC1~ADAD1', where 'SHC1' is the library gene, and 'ADAD1' is the anchor gene.

        Both genes will be targetFromSource AND interactingTargetFromSource, while keeping their roles.
        """
        genes = gene_pair.split('~')
        roles = [
            gene_role.replace('Combinations', '').lower(),
            'anchor'
        ]

        assert(len(genes) == 2)
        parsed = []

        for i, (gene, role) in enumerate(zip(genes, roles)):
            parsed.append({
                'targetFromSourceId': gene,
                'targetRole': role,
                'interactingTargetFromSourceId': genes[1] if i == 0 else genes[0],
                'interactingTargetRole': roles[1] if i == 0 else roles[0]
            })

        return parsed

    def get_lfc_data(self, lfc_file: str) -> DataFrame:
        """Reading log fold change file and generating a dataframe.

        The table has results from CRISPR/Cas9 experiments: the p-value, false discovery rate (fdr)
        and the log fold-change file describes the observed survival of the cells.

        Each row is a gene pair, columns contain results from the measurements. Column names contain
        information on the tested cell line (SIDM) and ocassionally the replicate (CSID) identifiers.

        This really wide table is stacked to get only three numeric columns (p-value, fold-change and FDR),
        plus string column with cell line.
        """

        # Fixed statistical field names:
        stats_fields = ['p-value', 'fdr', 'lfc']

        # Reading the data into a single dataframe:
        lfc_df = self.spark.read.csv(lfc_file, sep=' ', header=True)

        # Collect the cell lines from the lfc file header:
        cell_lines = set(['_'.join(x.split('_')[:-1]) for x in lfc_df.columns[4:]])

        # Generating struct for each cell lines:
        # SIDM00049_CSID1053_p-value
        # SIDM00049_CSID1053_fdr
        # SIDM00049_CSID1053_lfc
        # Into: SIDM00049_CSID1053: struct(p-value, fdr, lfc)
        expressions = map(
            lambda cell: (cell, struct([col(f'{cell}_{x}').cast(FloatType()).alias(x) for x in stats_fields])),
            cell_lines
        )

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, lfc_df)

        # Stack the previously generated columns:
        unpivot_expression = f'''stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)'''

        return (
            res_df

            # Unpivot:
            .select('id', 'Note1', 'Note2', expr(unpivot_expression))

            # Extracting and renaming log-fold-change parameters:
            .select(
                'id', 'cellLineName', 'Note1', 'Note2',
                col('cellLineData.lfc').alias('phenotypicConsequenceLogFoldChange'),
                col('cellLineData.p-value').alias('phenotypicConsequencePValue'),
                col('cellLineData.fdr').alias('phenotypicConsequenceFDR')
            )
        )

    def get_bliss_data(self, bliss_file: str) -> DataFrame:
        """Reading the BLISS cooperative effect file and generating a dataframe.

        This function is still being developed, as the format is likely subject to change.
        """
        # Fixed statistical field names:
        stats_fields = ['zscore', 'pval']

        # Read data:
        bliss_df = self.spark.read.csv(bliss_file, sep='\t', header=True)

        # Collect the cell lines from the lfc file header:
        cell_lines = set(['_'.join(x.split('_')[0:2]) for x in bliss_df.columns[4:] if x.startswith('SID')])

        # Generating struct for each cell lines:
        expressions = map(lambda cell: (cell, struct([col(f'{cell}_{x}').alias(x) for x in stats_fields])), cell_lines)

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, bliss_df)

        # Stack the previously generated columns:
        unpivot_expression = f'''stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)'''

        return (
            res_df

            # Create a consistent id column:
            .withColumn('id', regexp_replace(col('Gene_Pair'), ';', '~'))

            # Unpivot:
            .select('id', expr(unpivot_expression))

            # Extracting and renaming bliss statistical values:
            .select(
                'id', regexp_replace('cellLineName', '_strong', '').alias('cellLineName'),
                col('cellLineData.zscore').cast(FloatType()).alias('geneticInteractionScore'),
                col('cellLineData.pval').cast(FloatType()).alias('geneticInteractionPValue'),
            )
            .withColumn('statisticalMethod', lit('bliss'))
        )

    def get_gemini_data(self, gemini_file: str) -> DataFrame:
        """Parsing cooperative effects calculated by gemini method

        The structure of the file similar to the lfc file, the process follows a similar fashion,
        except it slightly different. Different enough to get a separate function.
        """

        # Fixed statistical field names:
        stats_fields = ['score', 'pval', 'FDR']

        # Reading the data into a single dataframe:
        gemini_df = self.spark.read.csv(gemini_file, sep=' ', header=True)

        # Collect the cell lines from the lfc file header:
        cell_lines = set(['_'.join(x.split('_')[:-1]) for x in gemini_df.columns[4:] if x.startswith('SID')])

        # There are some problems in joining gemini files on Encore side. It causes a serious issues:
        # 1. Multiple Gene_Pair columns in the file -> these will be indexed in the pyspark dataframe
        # 2. Some columns for some cell lines will be missing eg. pvalue for SIDM00049_CSID1053
        #
        # To mitigate these issue we have to check for gene pair header and remove cell lines with incomplete data.
        if 'Gene_Pair' in gemini_df.columns:
            gene_column = 'Gene_Pair'
        elif 'Gene_Pair0' in gemini_df.columns:
            gene_column = 'Gene_Pair0'
        else:
            raise ValueError(f'No Gene_Pair column in Gemini data: {",".join(gemini_df.columns)}')

        # We check if all stats columns available for all cell lines (this coming from a data joing bug at encore):
        missing_columns = [f'{cell}_{stat}' for cell in cell_lines for stat in stats_fields if f'{cell}_{stat}' not in gemini_df.columns]
        cells_to_drop = set(['_'.join(x.split('_')[:-1]) for x in missing_columns])

        # If there are missingness, the relevant cell lines needs to be removed from the analysis:
        if missing_columns:
            logging.warn(f'Missing columns: {", ".join(missing_columns)}')
            logging.warn(f'Dropping cell_lines: {", ".join(cells_to_drop)}')

            # Removing missing cell lines:
            cell_lines = [x for x in cell_lines if x not in cells_to_drop]

        # Generating struct for each cell lines:
        expressions = map(lambda cell: (cell, struct([col(f'{cell}_{x}').alias(x) for x in stats_fields])), cell_lines)

        # Applying map on the dataframe:
        res_df = reduce(lambda DF, value: DF.withColumn(*value), expressions, gemini_df)

        # Stack the previously generated columns:
        unpivot_expression = f'''stack({len(cell_lines)}, {", ".join([f"'{x}', {x}" for x in cell_lines])} ) as (cellLineName, cellLineData)'''

        return (
            res_df

            # Create a consistent id column:
            .withColumn('id', regexp_replace(col(gene_column), ';', '~'))

            # Unpivot:
            .select('id', expr(unpivot_expression))

            # Extracting and renaming gemini statistical values:
            .select(
                'id', regexp_replace('cellLineName', '_strong', '').alias('cellLineName'),
                col('cellLineData.score').cast(FloatType()).alias('geneticInteractionScore'),
                col('cellLineData.pval').cast(FloatType()).alias('geneticInteractionPValue'),
                col('cellLineData.FDR').cast(FloatType()).alias('geneticInteractionFDR')
            )
            .withColumn('statisticalMethod', lit('gemini'))
            .persist()
        )

    def parse_experiment(self, parameters: dict) -> DataFrame:
        """Parsing experiments based on the experimental descriptions.

        Args:
            parameters: Dictionary of experimental parameters. The following keys are required:
                - dataset: Name of the dataset eg. COLO1- referring to the first libraryset of colo.
                - diseaseFromSource: Name of the disease model of the experiment.
                - diseaseFromSourceMappedId: EFO ID of the disease model of the experiment.
                - logFoldChangeFile: File path to the log fold change file.
                - geminiFile: File path to the gemini file.
                - blissFile: File path to the bliss file.

        Returns:
            A pyspark dataframe of experiment data.

        Process:
            - reading all files.
            - Joining files.
            - Applying filters.
            - Adding additional columns + finalizing evidence model
        """

        disease_from_source = parameters['diseaseFromSource']
        disease_from_source_mapped_id = parameters['diseaseFromSourceMappedId']
        dataset = parameters['dataset']
        log_fold_change_file = parameters['logFoldChangeFile']
        gemini_file = parameters['geminiFile']
        # blissFile = parameters['blissFile'] # <= not used

        # Testing if experiment needs to be skipped:
        if parameters['skipStudy'] is True:
            logging.info(f'Skipping study: {dataset}')
            return None

        logging.info(f'Parsing experiment: {dataset}')

        # if no log fold change file is provided, we will not generate any evidence.
        if log_fold_change_file is None:
            logging.warning(f'No log fold change file provided for {dataset}.')
            return None

        # if no gemini file is provided, we will not generate any evidence.
        if gemini_file is None:
            logging.warning(f'No gemini file provided for {dataset}.')
            return None

        # Reading lfc data:
        lfc_file = f'{self.data_folder}/{log_fold_change_file}'
        lfc_df = self.get_lfc_data(lfc_file)
        logging.info(f'Number of gene paris in the log(fold change) dataset: {lfc_df.select("id").count()} rows')
        logging.info(f'Number cell lines in the log(fold change) dataset: {lfc_df.select("cellLineName").distinct().count()} rows')

        # Reading gemini data:
        gemini_file = f'{self.data_folder}/{gemini_file}'
        gemini_df = self.get_gemini_data(gemini_file)
        logging.info(f'Number of gene paris in the gemini dataset: {gemini_df.select("id").count()} rows')
        logging.info(f'Number cell lines in the gemini dataset: {gemini_df.select("cellLineName").distinct().count()} rows')

        """WARNING: Based on the communication with the encore team, the BLISS dataset is currently considered unreliable."""
        # bliss_file = f'{self.data_folder}/{blissFile}'
        # bliss_df = self.get_bliss_data(bliss_file)

        # Merging lfc + gemini:
        merged_dataset = (
            lfc_df

            # Data is joined by the gene-pair and cell line:
            # .join(bliss_df, how='inner', on=['id', 'cellLineName'])
            .join(gemini_df, how='inner', on=['id', 'cellLineName'])

            # Applying filters on logFoldChange + interaction p-value thresholds:
            .filter(
                (col('phenotypicConsequencePValue') <= self.log_fold_change_cutoff_p_val) &
                (col('phenotypicConsequenceFDR') <= self.log_fold_change_cutoff_FDR) &
                (col('geneticInteractionPValue') <= self.interaction_cutoff_p_val)
            )

            # Cleaning the cell line annotation:
            .withColumn('cellId', split(col('cellLineName'), '_').getItem(0))

            # Joining with cell passport data containing diseaseCellLines and biomarkers info:
            .join(
                self.cell_passport_df.select(col('id').alias('cellId'), 'diseaseCellLines', 'biomarkerList'),
                on='cellId', how='left'
            )
            .persist()
        )
        logging.info(f'Number of gene paris in the merged dataset: {merged_dataset.select("id").count()} rows')
        logging.info(f'Number cell lines in the merged dataset: {merged_dataset.select("cellLineName").distinct().count()} rows')

        evidence_df = (
            merged_dataset

            .withColumn('diseaseCellLines', array(col('diseaseCellLines')))

            # Parsing/exploding gene names and target roles:
            .withColumn('id', self.parse_targets(col('id'), col('Note1')))
            .select('*', explode(col('id')).alias('genes'))
            .select('*', col('genes.*'))

            # Adding some literals specific for this type of evidence:
            .withColumn('datatypeId', lit('ot_partner'))
            .withColumn('datasourceId', lit('encore'))
            .withColumn('projectId', lit('OTAR2062'))
            .withColumn('projectDescription', lit('Encore project'))  ### TODO - fix this!
            .withColumn('geneInteractionType', lit('cooperative'))

            # Adding disease information:
            .withColumn('diseaseFromSourceMappedId', lit(disease_from_source_mapped_id))
            .withColumn('diseaseFromSource', lit(disease_from_source))
            .withColumn('releaseVersion', lit(self.release_version))
            .withColumn('releaseDate', lit(self.release_date))

            # Removing unused columns:
            .drop(*['cellLineName', 'cellId', 'Note1', 'Note2', 'id', 'genes'])
            .distinct()

            .persist()
        )

        # If there's a warning message for the sudy, add that:
        if parameters['warningMessage'] is not None:
            evidence_df = evidence_df.withColumn('warningMessage', lit(parameters['warningMessage']))

        return evidence_df


def main(outputFile: str, parameters: dict, cellPassportFile: str) -> None:

    # Initialize spark session
    spark = initialize_sparksession()

    # Opening and parsing the cell passport data from Sanger:
    cell_passport_df = GenerateDiseaseCellLines(cellPassportFile, spark).get_mapping()

    logging.info(f'Cell passport dataframe has {cell_passport_df.count()} rows.')
    logging.info('Parsing experiment data...')

    # Initialising evidence generator:
    evidenceGenerator = EncoreEvidenceGenerator(spark, cell_passport_df, parameters['sharedMetadata'])

    # Create evidence for all experiments. Dataframes are collected in a list:
    evidence_dfs = []
    for experiment in parameters['experiments']:
        evidence_dfs.append(evidenceGenerator.parse_experiment(experiment))

    # Filter out None values, so only dataframes with evidence are kept:
    evidence_dfs = list(filter(None, evidence_dfs))

    # combine all evidence dataframes into one:
    combined_evidence_df = reduce(lambda df1, df2: df1.union(df2), evidence_dfs)

    # The combined evidence is written to a file:
    write_evidence_strings(combined_evidence_df, outputFile)


if __name__ == '__main__':

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(
        description='This script parses ENCORE data and generates disease-target evidence.')
    parser.add_argument('--output_file', '-o', type=str,
                        help='Output file. gzipped JSON', required=True)
    parser.add_argument('--parameter_file', '-p', type=str,
                        help='A JSON file describing exeriment metadata', required=True)
    parser.add_argument('--data_folder', '-d', type=str,
                        help='A folder in which the data files are located', required=True)
    parser.add_argument('--log_file', type=str,
                        help='File into which the logs are saved', required=False)
    parser.add_argument('--cell_passport_file', '-c', type=str,
                        help='File containing cell passport data', required=True)
    args = parser.parse_args()

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Validating data folder:
    if not os.path.isdir(args.data_folder):
        logging.error(f'Data folder {args.data_folder} does not exist.')
        raise ValueError(f'Data folder {args.data_folder} does not exist.')

    # Reading parameter json:
    try:
        with open(args.parameter_file, 'r') as f:
            parameters = json.load(f)
    except Exception as e:
        raise e(f'Could not read parameter file. {args.parameter_file}')

    # Updating parameters:
    parameters['sharedMetadata']['data_folder'] = args.data_folder

    # Passing all the required arguments:
    main(args.output_file, parameters, args.cell_passport_file)
