import logging
import os
import random
import time

from numpy import nan
from ontoma.interface import OnToma
from pandarallel import pandarallel
from pyspark.sql.functions import col, when
from pyspark.sql.types import StringType, StructField, StructType

ONTOMA_MAX_ATTEMPTS = 3
pandarallel.initialize()


def _simple_retry(func, **kwargs):
    """Simple retry handling for functions. Cannot be a decorator, so that the functions could still be pickled."""
    for attempt in range(1, ONTOMA_MAX_ATTEMPTS + 1):
        try:
            return func(**kwargs)
        except:
            # If this is not the last attempt, wait until the next one.
            if attempt != ONTOMA_MAX_ATTEMPTS:
                time.sleep(5 + 10 * random.random())
    logging.error(f'OnToma lookup failed for {kwargs!r}')
    return []


def _ontoma_udf(row, ontoma_instance):
    """Try to map first by disease name (because that branch of OnToma is more stable), then by disease ID."""
    disease_name = None
    if row['diseaseFromSource']:
        disease_name = ' '.join(row['diseaseFromSource'].replace('obsolete', '').split())
    disease_id = row['diseaseFromSourceId'].replace('_', ':') if row['diseaseFromSourceId'] else None
    mappings = []
    if disease_name:
        mappings = _simple_retry(ontoma_instance.find_term, query=disease_name, code=False)
    if not mappings and disease_id and ':' in disease_id:
        mappings = _simple_retry(ontoma_instance.find_term, query=disease_id, code=True)
    return [m.id_ot_schema for m in mappings]


def add_efo_mapping(evidence_strings, spark_instance, ontoma_cache_dir=None, efo_version=None):
    """Given evidence strings with diseaseFromSource and diseaseFromSourceId fields, try to populate EFO mapping
    field diseaseFromSourceMappedId. In case there are multiple matches, the evidence strings will be exploded
    accordingly.

    Currently, both source columns (diseaseFromSource and diseaseFromSourceId) need to be present in the original
    schema, although they do not have to be populated for all rows."""
    logging.info('Collect all distinct (disease name, disease ID) pairs.')
    disease_info_to_map = evidence_strings.select('diseaseFromSource', 'diseaseFromSourceId').distinct().toPandas()

    # If no EFO version is specified:
    if not efo_version:
        # try to extract from environment variable.
        if "EFO_VERSION" in os.environ:
            efo_version = os.environ["EFO_VERSION"]
        # Set default version to latest.
        else:
            logging.warning('No EFO version specified. Using latest version.')
            efo_version = 'latest'

    logging.info(f'Initialise OnToma instance. Using EFO version {efo_version}')
    ontoma_instance = OnToma(cache_dir=ontoma_cache_dir, efo_release=efo_version)

    logging.info('Map disease information to EFO.')
    disease_info_to_map['diseaseFromSourceMappedId'] = disease_info_to_map.parallel_apply(
        _ontoma_udf, args=(ontoma_instance,), axis=1
    )
    disease_info_to_map = (
        disease_info_to_map.explode('diseaseFromSourceMappedId')
        # Cast all null values to python None to avoid errors in Spark's DF
        .fillna(nan).replace([nan], [None])
    )

    logging.info('Join the resulting information into the evidence strings.')
    schema = StructType(
        [
            StructField("diseaseFromSource_right", StringType(), True),
            StructField("diseaseFromSourceId_right", StringType(), True),
            StructField("diseaseFromSourceMappedId", StringType(), True),
        ]
    )
    disease_info_df = spark_instance.createDataFrame(disease_info_to_map, schema=schema).withColumn(
        'diseaseFromSourceMappedId', when(col('diseaseFromSourceMappedId') != 'nan', col('diseaseFromSourceMappedId'))
    )
    # WARNING: Spark's join operator is not null safe by default and most of the times, `diseaseFromSourceId` will be null.
    # `eqNullSafe` is a special null safe equality operator that is used to join the two dataframes.
    join_cond = (evidence_strings.diseaseFromSource.eqNullSafe(disease_info_df.diseaseFromSource_right)) & (
        evidence_strings.diseaseFromSourceId.eqNullSafe(disease_info_df.diseaseFromSourceId_right)
    )
    return evidence_strings.join(disease_info_df, on=join_cond, how='left').drop(
        'diseaseFromSource_right', 'diseaseFromSourceId_right'
    )
