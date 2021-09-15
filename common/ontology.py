import logging
import random
import time

from ontoma.interface import OnToma
from pandarallel import pandarallel
from pyspark.sql.functions import col, when

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
    disease_name = row['diseaseFromSource']
    disease_id = row['diseaseFromSourceId'].replace('_', ':') if row['diseaseFromSourceId'] else None
    mappings = []
    if disease_name:
        mappings = _simple_retry(ontoma_instance.find_term, query=disease_name, code=False)
    if not mappings and disease_id and ':' in disease_id:
        mappings = _simple_retry(ontoma_instance.find_term, query=disease_id, code=True)
    return [m.id_ot_schema for m in mappings]


def add_efo_mapping(evidence_strings, spark_instance, ontoma_cache_dir=None):
    """Given evidence strings with diseaseFromSource and diseaseFromSourceId fields, try to populate EFO mapping
    field diseaseFromSourceMappedId. In case there are multiple matches, the evidence strings will be exploded
    accordingly.

    Currently, both source columns (diseaseFromSource and diseaseFromSourceId) need to be present in the original
    schema, although they do not have to be populated for all rows."""
    logging.info('Collect all distinct (disease name, disease ID) pairs.')
    disease_info_to_map = (
        evidence_strings
        .select('diseaseFromSource', 'diseaseFromSourceId')
        .distinct()
        .toPandas()
    )

    logging.info('Initialise OnToma instance')
    ontoma_instance = OnToma(cache_dir=ontoma_cache_dir)

    logging.info('Map disease information to EFO.')
    disease_info_to_map['diseaseFromSourceMappedId'] = disease_info_to_map.parallel_apply(
        _ontoma_udf, args=(ontoma_instance,), axis=1
    )
    disease_info_to_map = disease_info_to_map.explode('diseaseFromSourceMappedId')

    logging.info('Join the resulting information into the evidence strings.')
    disease_info_df = (
        spark_instance
        .createDataFrame(disease_info_to_map.astype(str))
        .withColumn(
            'diseaseFromSourceMappedId',
            when(col('diseaseFromSourceMappedId') != 'nan', col('diseaseFromSourceMappedId'))
        )
    )
    return evidence_strings.join(
        disease_info_df,
        on=['diseaseFromSource', 'diseaseFromSourceId'],
        how='left'
    )
