import argparse
import logging

from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F

from common.evidence import convert_stringified_array_to_array, load_sql_table_to_spark, write_evidence_strings

# Define data of interest
PROBES_SETS = [
    'Gray Laboratory Probes',
    'Bromodomains chemical toolbox',
    'Chemical Probes.org',
    'Open Science Probes',
    'High-quality chemical probes',
    'opnMe Portal',
    'High-quality chemical probes',
    'Probe Miner (suitable probes)',
    'Protein methyltransferases chemical toolbox',
    'SGC Probes',
    'Natural product-based probes and drugs',
    'JUMP-Target 1 Compound Set',
    'JUMP-Target 2 Compound Set',
    'Chemical Probes for Understudied Kinases',
]
SCORES = [
    'Probe Miner Score',
    'Cells score (Chemical Probes.org)',
    'Organisms score (Chemical Probes.org)',
    'P&D probe-likeness score',
]
ORGANISMS = ['Homo sapiens (Human)', 'Homo sapiens']

# DB params
PATH_TO_JAR_FILE = '/Users/irene/MEGAsync/EBI/repos/evidence_datasource_parsers/postgresql-42.4.0.jar'
DB_URL = 'jdbc:postgresql://localhost:5432/probes'
DB_USER = 'irene'
DB_PASSWORD = 'probes'


def main(db_url: str, db_user: str, db_password: str, drug_index: str):

    # Load data
    compound_df = (
        load_sql_table_to_spark(db_url, db_user, db_password, 'compound', ['compoundid', 'name', 'inchikey'])
        .withColumnRenamed('compoundid', 'compound_id')
        .withColumnRenamed('name', 'compound_name')
    )
    probe_df = (
        load_sql_table_to_spark(
            db_url, db_user, db_password, 'probe', ['probeid', 'compound_id', 'control', 'origin_id', 'obsolete_flag']
        )
        .withColumnRenamed('probeid', 'probe_id')
        .filter(F.col('obsolete_flag') == 't')
        .drop('obsolete_flag')
    )
    target_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'target', ['targetid', 'organism_id', 'uniprotid']
    ).withColumnRenamed('targetid', 'target_id')
    compound_set_df = (
        load_sql_table_to_spark(db_url, db_user, db_password, 'compoundset', ['compoundsetid', 'name', 'source_url'])
        .withColumnRenamed('compoundsetid', 'compoundset_id')
        .withColumnRenamed('name', 'source_name')
        .filter(F.col('name').isin(PROBES_SETS))
    )
    compound_to_compoundset_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'compoundtocompoundset', ['compound_id', 'compoundset_id']
    )
    probe_to_base_target_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'probetobasetarget', ['probe_id', 'basetarget_id']
    )
    base_target_to_target_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'targettobasetarget', ['basetarget_id', 'target_id']
    )
    score_df = (
        load_sql_table_to_spark(db_url, db_user, db_password, 'score', ['scoreid', 'name'])
        .select(F.col('scoreid').alias('score_id'), F.col('name').alias('score_name'))
        .filter(F.col('score_name').isin(SCORES))
    )
    probe_target_score_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'compoundtargetscore', ['compound_id', 'basetarget_id', 'percentage', 'score_id']
    ).withColumnRenamed('percentage', 'score_percentage')
    compound_action_df = load_sql_table_to_spark(
        db_url, db_user, db_password, 'compoundaction', ['compound_id', 'target_id', 'actiontype_id']
    )
    organism_df = (
        load_sql_table_to_spark(db_url, db_user, db_password, 'organism', ['organismid', 'name'])
        .withColumnRenamed('organismid', 'organism_id')
        .filter(F.col('name').isin(ORGANISMS))
        .drop('name')
    )
    drug_df = spark.read.parquet(drug_index).select(F.col('id').alias('drugId'), F.col('inchiKey').alias('inchikey'))

    # Join data
    cols_to_keep = [
        'uniprotid',
        'compound_name',
        'control',
        'origin_id',
        'drugId',
        'inchikey',
        'actiontype_id',
        'score_name',
        'score_percentage',
        'source_name',
        'source_url',
    ]
    probe_df = (
        # Linkage between the probe and the target
        probe_df.join(probe_to_base_target_df, on='probe_id', how='inner')
        .join(base_target_to_target_df, on='basetarget_id', how='inner')
        .join(target_df, on='target_id', how='inner')
        .join(organism_df, on='organism_id', how='inner')
        # Filter probes on compound sets of interest
        .join(compound_to_compoundset_df, on='compound_id', how='inner')
        .join(compound_set_df, on='compoundset_id', how='inner')
        # Extract compound information
        .join(compound_df, on='compound_id', how='inner')
        .join(drug_df, on='inchikey', how='left')
        .join(compound_action_df, on=['compound_id', 'target_id'], how='left')
        # Extract score information
        .join(probe_target_score_df, on=['compound_id', 'basetarget_id'], how='left')
        .join(score_df, on='score_id', how='left')
        .filter(F.col('score_name').isin(SCORES))
        .select(cols_to_keep)
        .distinct()
        .persist()
    )

    probe_df = handle_high_quality_probes(probe_df)

    return transform_data(probe_df)


def handle_high_quality_probes(probe_df: DataFrame):
    '''
    One of the sources collects probes of high quality based on a series of parameters defined by P&Ds.
    We want to report these not as a separate probes set, but as a flag in the probe table that will be present in all sources.
    For that:
        - `isHighQuality` is set to `true` if the probe is in the HQ list
        - the records in the HQ set will be removed not to contain duplicate information
    There's the caveat of the fact that some probes are only reported in the HQ source. We will keep these and will be reported as P&Ds as a source.
    '''

    hq_probes = (
        probe_df.filter(F.col('source_name').isin('High-quality chemical probes'))
        .select('compound_name')
        .distinct()
        .toPandas()['compound_name']
        .to_list()
    )
    only_hq_probes = (
        probe_df.filter(F.col('compound_name').isin(hq_probes))
        .groupby('compound_name')
        .agg(F.collect_set('source_name').alias('source_names'))
        .filter(F.size('source_names') == 1)
        .filter(F.array_contains('source_names', 'High-quality chemical probes'))
        .select('compound_name')
        .distinct()
        .toPandas()['compound_name']
        .to_list()
    )
    logging.info(f'{len(only_hq_probes)} probes are only reported in the HQ source.')
    logging.info(f'{len(hq_probes)} probes are flagged as high quality.')

    return (
        probe_df.withColumn('isHighQuality', F.when(F.col('compound_name').isin(hq_probes), True).otherwise(False))
        .withColumn(
            'source_name',
            F.when(F.col('compound_name').isin(only_hq_probes), 'High-quality chemical probes').otherwise(
                F.col('source_name')
            ),
        )
        .filter(F.col('source_name') != 'High-quality chemical probes')
    )


def transform_data(probe_df: DataFrame):
    '''
    Applies transformations to the data to follow the chemical probes data model.
    '''

    probe_df = (
        probe_df
        # Clean column that collects the probes controls
        .transform(lambda df: convert_stringified_array_to_array(df, 'control'))
        .withColumn('control', F.explode_outer('control'))
        .filter(F.col('control').isNull() | ~F.col('control').isin(['1', 'No control']))
        # Build a more complete URL
        .withColumn(
            'source_url',
            F.when(F.col('source_url').contains('probeminer'), F.concat('source_url', 'uniprotid'))
            .when(
                F.col('source_url').contains('probes.org'),
                F.concat(F.lit('https://www.chemicalprobes.org/?q='), 'compound_name'),
            )
            .when(F.col('source_url').contains('thesgc'), F.concat('source_url', F.lit('/'), 'compound_name'))
            .when(F.col('source_url').contains('graylab'), F.lit(None))
            .when(
                F.col('source_url').contains('opnme'),
                F.concat(F.lit('https://opnme.com/search-results/'), 'compound_name'),
            )
            .when(
                F.col('source_url').contains('frankfurt'),
                F.concat(
                    F.element_at(F.split(F.col('source_url'), '#!start'), 1),
                    F.lit('#!specificprobeoverview/'),
                    'compound_name',
                ),
            )
            .otherwise(F.col('source_url')),
        )
        # Collect the source of the probe in a struct
        .withColumn('url', F.struct(F.col('source_name').alias('niceName'), F.col('source_url').alias('url')))
        # Split scores by type
        .withColumn('probeMinerScore', F.when(F.col('score_name') == 'Probe Miner Score', F.col('score_percentage')))
        .withColumn(
            'scoreInCells',
            F.when(F.col('score_name') == 'Cells score (Chemical Probes.org)', F.col('score_percentage')),
        )
        .withColumn(
            'scoreInOrganisms',
            F.when(F.col('score_name') == 'Organisms score (Chemical Probes.org)', F.col('score_percentage')),
        )
        .withColumn(
            'probesDrugsScore', F.when(F.col('score_name') == 'P&D probe-likeness score', F.col('score_percentage'))
        )
        .groupby('compound_name', 'uniprotid', 'control', 'drugId')
        .agg(
            F.collect_set('origin_id').alias('origin'),
            F.collect_set('isHighQuality').alias('isHighQuality'),
            F.collect_set('actiontype_id').alias('mechanismOfAction'),
            F.collect_set('probeMinerScore').alias('probeMinerScore'),
            F.collect_set('scoreInCells').alias('scoreInCells'),
            F.collect_set('scoreInOrganisms').alias('scoreInOrganisms'),
            F.collect_set('probesDrugsScore').alias('probesDrugsScore'),
            F.collect_set('url').alias('urls'),
        )
        # The max score is extracted in those cases that there is more than one value - principally ProbeMiner
        .withColumn('probeMinerScore', F.array_max('probeMinerScore'))
        .withColumn('scoreInCells', F.array_max('scoreInCells'))
        .withColumn('scoreInOrganisms', F.array_max('scoreInOrganisms'))
        .withColumn('probesDrugsScore', F.array_max('probesDrugsScore'))
        # Convert empty arrays into null
        .withColumn(
            'mechanismOfAction',
            F.when(F.size(F.col('mechanismOfAction')) == 0, F.lit(None)).otherwise(F.col('mechanismOfAction')),
        )
        .withColumn('origin', F.when(F.size(F.col('origin')) == 0, F.lit(None)).otherwise(F.col('origin')))
        .withColumn('isHighQuality', F.when(F.array_contains(F.col('isHighQuality'), True), True).otherwise(False))
        .withColumnRenamed('uniprotid', 'targetFromSourceId')
        .withColumnRenamed('compound_name', 'id')
        .distinct()
    )

    logging.info(f'{probe_df.count()}  have been processed.')

    return probe_df


def get_parser():
    'Get parser object for script ChemicalProbes.py.'
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '--probes_db_path',
        help='Path to the Probes&Drugs PostgreSQL db to establish a connection to read the data on Spark.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--drug_index',
        help='Directory of parquet files that stores OT\'s disease index.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Output gzipped json file following the chemical probes data model.',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--log_file',
        help='Destination of the logs generated by this script. Defaults to None',
        type=str,
        nargs='?',
        default=None,
    )

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    # Logger initializer. If no log_file is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    spark = (
        SparkSession.builder.appName('app')
        .config("spark.jars", PATH_TO_JAR_FILE)
        .config('spark.driver.extraClassPath', PATH_TO_JAR_FILE)
        .getOrCreate()
    )

    probe_df = main(db_url=args.probes_db_path, db_user=DB_USER, db_password=DB_PASSWORD, drug_index=args.drug_index)
    write_evidence_strings(probe_df, args.output)
    logging.info(f'Probes dataset has been saved to {args.output}. Exiting.')
