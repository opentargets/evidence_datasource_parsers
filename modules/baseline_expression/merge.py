#!/usr/bin/env python3
"""
Merges multiple parquet datasets into a single parquet output and (optionally) writes a small JSON sample.

Requires Apache Spark (pyspark) for distributed processing.

Expected input layout (per dataset):
- Base directory: <base>/<aggregation>/<dataset>/parquet/*
  where <aggregation> is one of: aggregated, unaggregated
  and <dataset> is a directory name (e.g., gtex, pride, tabula_sapiens, dice)

Outputs:
- Merged parquet: <base>/<aggregation>/merged/parquet
- Small validation sample (JSON, optional, disable with --no-sample): <base>/<aggregation>/merged_sample/json
"""

import argparse
import os

from pyspark.sql import SparkSession
from pyspark.sql.functions import col, rand, row_number
from pyspark.sql.window import Window


class MergeParquetDatasets:
    """Collection of steps to read, merge, repartition, and write parquet datasets."""

    def read_input_data(self):
        """Read the input parquet datasets and union them with schema alignment."""
        paths = [
            f"{self.base_directory_path.rstrip('/')}/{self.aggregation}/{d}/parquet/*"
            for d in self.datasets
        ]

        if not paths:
            raise ValueError("No dataset paths were constructed. Check --datasets.")

        print("Reading parquet inputs:")
        for p in paths:
            print(f"  • {p}")

        # Read first dataset
        df = self.spark.read.parquet(paths[0])

        # Union remaining datasets, aligning by column name
        for p in paths[1:]:
            df_next = self.spark.read.parquet(p)
            df = df.unionByName(df_next, allowMissingColumns=True)

        self.df = df
        print("Finished reading and merging input data.")

    def pack_data_for_output(self):
        """Repartition and write merged parquet; optionally write a small JSON sample."""
        # Resolve output base (optionally prefixed with file:// for local mode)
        def _out(path):
            if self.local:
                return f"file://{path}"
            return path

        merged_parquet_out = (
            self.output_directory_path
            if self.output_directory_path
            else f"{self.base_directory_path.rstrip('/')}/{self.aggregation}/merged/parquet"
        )
        sample_json_out = (
            self.sample_output_directory_path
            if self.sample_output_directory_path
            else f"{self.base_directory_path.rstrip('/')}/{self.aggregation}/merged_sample/json"
        )

        # Count rows and decide number of output files
        print("Counting rows...")
        total_rows = self.df.count()
        print(f"Total rows detected: {total_rows:,}")

        if self.num_output_files is not None and self.num_output_files > 0:
            num_files = int(self.num_output_files)
        else:
            rpf = int(self.rows_per_file) if self.rows_per_file and self.rows_per_file > 0 else 10_000_000
            # ceil division, minimum 1
            num_files = max(1, (total_rows + rpf - 1) // rpf)

        print(f"Target number of output files: {num_files}")

        # Repartition for balanced file sizes
        print("Repartitioning (shuffle)...")
        df_out = self.df.repartition(num_files)

        # Write merged parquet
        print(f"Writing merged parquet to: {merged_parquet_out}")
        df_out.write.mode("overwrite").parquet(_out(merged_parquet_out))
        print("Finished writing merged parquet output.")

        # Optionally write small JSON sample, balanced across datasourceId if present
        if self.no_sample:
            print("Skipping sample JSON (--no-sample set).")
            return

        n_rows = int(self.sample_rows) if self.sample_rows and self.sample_rows > 0 else 0
        if n_rows <= 0:
            print("Skipping sample JSON (sample_rows <= 0).")
            return

        print(f"Preparing JSON sample of ~{n_rows} rows...")
        cols = set(self.df.columns)
        if "datasourceId" in cols:
            n_datasources = self.df.select("datasourceId").distinct().count()
            # ceil division
            k = (n_rows + n_datasources - 1) // max(1, n_datasources)
            w = Window.partitionBy("datasourceId").orderBy(rand())
            sample_df = (
                self.df.withColumn("rn", row_number().over(w))
                       .where(col("rn") <= k)
                       .drop("rn")
                       .limit(n_rows)
            )
        else:
            # Fallback: simple random sample if no datasourceId column
            print("Column 'datasourceId' not found; sampling without stratification.")
            sample_df = self.df.orderBy(rand()).limit(n_rows)

        print(f"Writing JSON sample to: {sample_json_out}")
        sample_df.write.mode("overwrite").json(_out(sample_json_out))
        print("Finished writing sample JSON output.")

    def main(self):
        # Build Spark 
        self.spark = (
            SparkSession.builder
            .appName("MergeParquetDatasets")
            .config("spark.jars", "https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar")
            .config("spark.executor.memory", "70g")
            .config("spark.driver.memory", "50g")
            .config("spark.memory.offHeap.enabled", True)
            .config("spark.memory.offHeap.size", "16g")
            .config("spark.sql.files.ignoreCorruptFiles", "true")
            .config("spark.sql.files.ignoreMissingFiles", "true")
            .getOrCreate()
        )
        # Avoid 'HIVE_DEFAULT_PARTITION' for null partitions in some environments
        self.spark.conf.set("hive.exec.default.partition.name", "None")

        print("Starting merge pipeline...")
        self.read_input_data()
        print("Packing data for output...")
        self.pack_data_for_output()
        self.spark.stop()
        print("All done.")

    def __init__(
        self,
        base_directory_path: str,
        aggregation: str,
        datasets,
        output_directory_path: str | None = None,
        sample_output_directory_path: str | None = None,
        rows_per_file: int = 10_000_000,
        num_output_files: int | None = None,
        sample_rows: int = 100,
        local: bool = False,
        no_sample: bool = False,
    ):
        self.base_directory_path = base_directory_path
        self.aggregation = aggregation
        # Support comma-separated string or list
        if isinstance(datasets, str):
            self.datasets = [d.strip() for d in datasets.split(",") if d.strip()]
        else:
            self.datasets = list(datasets)

        self.output_directory_path = output_directory_path
        self.sample_output_directory_path = sample_output_directory_path
        self.rows_per_file = rows_per_file
        self.num_output_files = num_output_files
        self.sample_rows = sample_rows
        self.local = local
        self.no_sample = no_sample
        self.spark = None  # Will be initialized in main()


parser = argparse.ArgumentParser(description="Merge multiple parquet datasets and optionally emit a small JSON sample.")
parser.add_argument(
    "--local", action="store_true", default=False,
    help="Prefix outputs with file:// for local filesystem writes."
)
parser.add_argument(
    "--base-directory-path", required=False, type=str, default="/home/alegbe/results",
    help="Root directory containing <aggregation>/<dataset>/parquet/*"
)
parser.add_argument(
    "--aggregation", required=False, type=str, default="aggregated",
    choices=["aggregated", "unaggregated"],
    help="Which aggregation subdirectory to read from."
)
parser.add_argument(
    "--datasets", required=False, type=str, default="gtex,pride,tabula_sapiens,dice",
    help="Comma-separated list of dataset directory names to merge."
)
parser.add_argument(
    "--output-directory-path", required=False, type=str, default=None,
    help="Directory for merged parquet output. Defaults to <base>/<aggregation>/merged/parquet"
)
parser.add_argument(
    "--sample-output-directory-path", required=False, type=str, default=None,
    help="Directory for JSON sample output. Defaults to <base>/<aggregation>/merged_sample/json"
)
parser.add_argument(
    "--rows-per-file", required=False, type=float, default=1e7,
    help="Approximate rows per output file (ignored if --num-output-files is set)."
)
parser.add_argument(
    "--num-output-files", required=False, type=int, default=None,
    help="Force an exact number of output files (overrides --rows-per-file)."
)
parser.add_argument(
    "--sample-rows", required=False, type=int, default=100,
    help="Number of rows to include in the JSON sample."
)
parser.add_argument(
    "--no-sample", action="store_true", default=False,
    help="Disable writing the JSON sample (overrides --sample-rows)."
)

if __name__ == "__main__":
    args = parser.parse_args()
    MergeParquetDatasets(
        args.base_directory_path,
        args.aggregation,
        args.datasets,
        args.output_directory_path,
        args.sample_output_directory_path,
        rows_per_file=int(args.rows_per_file) if args.rows_per_file else 10_000_000,
        num_output_files=args.num_output_files,
        sample_rows=args.sample_rows,
        local=args.local,
        no_sample=args.no_sample,
    ).main()
