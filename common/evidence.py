from decimal import Decimal
import os
import tempfile

from psutil import virtual_memory
from pyspark.conf import SparkConf
from pyspark.sql import SparkSession


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


def get_exponent(number: float) -> int:
    """Get the exponent of a number."""
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1


def get_mantissa(number: float) -> float:
    """Get the mantissa of a number."""
    return round(number / (10 ** get_exponent(number)), 3)
