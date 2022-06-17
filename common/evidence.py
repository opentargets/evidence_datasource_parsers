import os
import tempfile

from psutil import virtual_memory
from pyspark import SparkFiles
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import col


def detect_spark_memory_limit():
    """Spark does not automatically use all available memory on a machine. When working on large datasets, this may
    cause Java heap space errors, even though there is plenty of RAM available. To fix this, we detect the total amount
    of physical memory and allow Spark to use (almost) all of it."""
    mem_gib = virtual_memory().total >> 30
    return int(mem_gib * 0.9)


def write_evidence_strings(evidence, output_file):
    """Exports the table to a compressed JSON file containing the evidence strings."""
    with tempfile.TemporaryDirectory() as tmp_dir_name:
        (
            evidence.coalesce(1)
            .write.format('json')
            .mode('overwrite')
            .option('compression', 'org.apache.hadoop.io.compress.GzipCodec')
            .save(tmp_dir_name)
        )
        json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
        assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
        os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)


def initialize_sparksession() -> SparkSession:
    """Initialize spark session."""

    spark_mem_limit = detect_spark_memory_limit()
    spark_conf = (
        SparkConf()
        .set('spark.driver.memory', f'{spark_mem_limit}g')
        .set('spark.executor.memory', f'{spark_mem_limit}g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        .set('spark.ui.showConsoleProgress', 'false')
    )
    spark = (
        SparkSession.builder.config(conf=spark_conf)
        .master('local[*]')
        .config("spark.driver.bindAddress", "127.0.0.1")
        .getOrCreate()
    )

    return spark


def read_path(path: str, spark_instance) -> DataFrame:
    """Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single TSV file; a single JSON file; a directory with JSON files; a directory with Parquet files."""
    if path is None:
        return None

    # The provided path must exist and must be either a file or a directory.
    assert os.path.exists(path), f'The provided path {path} does not exist.'
    assert os.path.isdir(path) or os.path.isfile(path), f'The provided path {path} is neither a file or a directory.'

    # Case 1: We are provided with a single file.
    if os.path.isfile(path):
        if path.endswith('.csv'):
            return spark_instance.read.csv(path, header=True, inferSchema=True)
        if path.endswith('.tsv'):
            return spark_instance.read.csv(path, sep='\t', header=True)
        elif path.endswith(('.json', '.json.gz', '.jsonl', '.jsonl.gz')):
            return spark_instance.read.json(path)
        else:
            raise AssertionError(f'The format of the provided file {path} is not supported.')

    # Case 2: We are provided with a directory. Let's peek inside to see what it contains.
    all_files = [os.path.join(dp, filename) for dp, dn, filenames in os.walk(path) for filename in filenames]

    # It must be either exclusively JSON, or exclusively Parquet.
    json_files = [fn for fn in all_files if fn.endswith(('.json', '.json.gz', '.jsonl', '.jsonl.gz'))]
    parquet_files = [fn for fn in all_files if fn.endswith('.parquet')]
    assert not (json_files and parquet_files), f'The provided directory {path} contains a mix of JSON and Parquet.'
    assert json_files or parquet_files, f'The provided directory {path} contains neither JSON nor Parquet.'

    # A directory with JSON files.
    if json_files:
        return spark_instance.read.option('recursiveFileLookup', 'true').json(path)

    # A directory with Parquet files.
    if parquet_files:
        return spark_instance.read.parquet(path)


def import_trait_mappings() -> DataFrame:
    """Load the remote trait mappings file to a Spark dataframe."""

    remote_trait_mappings_url = (
        'https://raw.githubusercontent.com/opentargets/curation/master/mappings/disease/manual_string.tsv'
    )

    SparkSession.getActiveSession().sparkContext.addFile(remote_trait_mappings_url)

    return (
        SparkSession.getActiveSession()
        .read.csv(SparkFiles.get('manual_string.tsv'), header=True, sep='\t')
        .select(
            col('PROPERTY_VALUE').alias('diseaseFromSource'), col('SEMANTIC_TAG').alias('diseaseFromSourceMappedId')
        )
    )
