"""Ingests GTEx V10 data and generates the unaggregated baseline expression data.

Requires Apache Spark (pyspark) for distributed processing.

Expected input formats:
- GTEx source data: gzipped GCT file (.gct.gz) containing sample-by-gene TPM counts.
- Sample metadata: TSV file with columns SAMPID, SMTSD, SMUBRID.
- Subject metadata: TSV file with columns SUBJID, AGE, SEX.
"""

import argparse
import functools
import gzip
import json
import os

from pyspark.sql.functions import split, col, array, struct, lit, explode, concat_ws, when
from pyspark.sql import SparkSession


class BaselineExpression:
    """Collection of steps to generate unaggregated baseline expression data."""

    def read_gtex_data(
        self,
    ):

        # Read the .gct.gz via RDD, skip first two lines
        rdd = self.spark.sparkContext.textFile(self.gtex_source_data_path)
        data_rdd = (
            rdd
            .zipWithIndex()
            .filter(lambda x: x[1] >= 2)  # drop lines 0 & 1
            .map(lambda x: x[0])
        )

        # Parse as TSV with header starting from third line
        df = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(data_rdd)
        )
        
        # Clean up columns
        df = (
            df
            .drop("Description")                      # drop unused col
            .filter(~col("Name").endswith("_PAR_Y"))  # remove genes ending in _PAR_Y
            .withColumn("Name", split(col("Name"), r"\.").getItem(0))  # strip .2 etc.
        )

        # If matrix mode, keep in wide format and store separately
        if self.matrix:
            # Store the wide format matrix
            data_cols = [c for c in df.columns if c != "Name"]
            self.df_matrix = df.withColumnRenamed("Name", "targetId")
            self.sample_ids = data_cols  # Keep track of sample column names
            # Don't convert to long format, just return
            return

        # Wide→long: explode an array of structs (one per sample column)
        data_cols = [c for c in df.columns if c != "Name"]

        # build an array<struct<Sample:string,TPM:double>>  
        samp_structs = array(*[
            struct(lit(c).alias("Sample"), col(c).cast("double").alias("TPM"))
            for c in data_cols
        ]).alias("samp_tpm")

        df_long = (
            df
            .select("Name", explode(samp_structs).alias("x"))
            .select(
                col("Name"),
                col("x.Sample").alias("OrigSample"),
                col("x.TPM")
            )
        )

        # split OrigSample into donorId + sampleId (split on second “-”)
        parts = split(col("OrigSample"), "-", 3)
        df_long = (
            df_long
            .withColumn("donorId", concat_ws("-", parts.getItem(0), parts.getItem(1)))
            .withColumn("sampleId", parts.getItem(2))
        )

        # Read sample metadata and rename columns
        meta_sample = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(self.sample_metadata_path)
            .select(
                col("SAMPID").alias("OrigSample"),
                col("SMTSD").alias("Tissue"),
                col("SMUBRID").alias("TissueOntologyID"),
            )
        )
        # Read subject metadata and rename columns, also map sex from 1/2 to M/F
        meta_subject = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(self.subject_metadata_path)
            .select(
                col("SUBJID").alias("donorId"),
                col("AGE").alias("Age"),
                col("SEX").alias("Sex")
            )
            .withColumn(
                "Sex",
                when(col("Sex") == "1", lit("M"))
                .when(col("Sex") == "2", lit("F"))
                .otherwise(lit("U"))
            )
        )

        # Join the long DataFrame with metadata
        df_long = (
            df_long
            .join(meta_sample,  on="OrigSample", how="left")
            .join(meta_subject, on="donorId",    how="left")
            .select(
                "Name", "donorId", "sampleId", "TPM",
                "Tissue", "TissueOntologyID", "Age", "Sex"
            )
            # Replace : with _ in TissueOntologyID
            .withColumn(
                "TissueOntologyID",
                when(col("TissueOntologyID").isNotNull(), 
                     concat_ws("_", split(col("TissueOntologyID"), ":")))
                .otherwise(lit(None))
            )
        )

        # # Renamed most columns to match the desired output format
        self.df = (
            df_long
            .withColumnRenamed("Name", "targetId")
            .withColumnRenamed("Tissue", "tissueBiosampleFromSource")
            .withColumnRenamed("TissueOntologyID", "tissueBiosampleId")
            .withColumn("unit", lit("TPM"))  # Add unit column with constant value "TPM"
            .withColumn("datasourceId", lit("gtex"))  # Add datasourceId column with constant value "gtex"
            .withColumn("datatypeId", lit("bulk rna-seq"))  # Add datatypeId column with constant value "bulk rna-seq"
            .withColumnRenamed("TPM", "expression")
            .withColumn("donorId", col("DonorID"))
            .withColumnRenamed("Sex", "sex")
            .withColumnRenamed("Age", "age")
        ).drop("sampleId")

    def pack_data_for_output(self, local: bool = False, json: bool = False, matrix: bool = False):
        """Use spark to write the DataFrame to parquet format."""
        if local:
            output_path = f"file://{self.output_directory_path}/"
        else:
            output_path = f"{self.output_directory_path}/"
        
        if matrix:
            # Save in matrix form with separate metadata
            self.save_as_matrix(output_path)
        elif json:
            output_path = f"{output_path}/json/gtex_tissue"
            # If JSON output is requested, convert DataFrame to JSON format
            self.df.write.mode("overwrite").json(output_path)
            print(f"Data written to {output_path} in JSON format")
        else:
            output_path = f"{output_path}/parquet/gtex_tissue"
            # If parquet output is requested, convert DataFrame to parquet format
            self.df.write.mode("overwrite").parquet(output_path)
            print(f"Data written to {output_path}")

    def save_as_matrix(self, output_path: str):
        """Save expression data in matrix form with separate metadata file."""
        import pandas as pd
        
        # Matrix is already in wide format from read_gtex_data
        # Just need to prepare metadata from the sample IDs
        
        # Read sample metadata
        meta_sample = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(self.sample_metadata_path)
            .select(
                col("SAMPID").alias("OrigSample"),
                col("SMTSD").alias("tissueBiosampleFromSource"),
                col("SMUBRID").alias("TissueOntologyID"),
            )
        )
        
        # Read subject metadata and map sex from 1/2 to M/F
        meta_subject = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(self.subject_metadata_path)
            .select(
                col("SUBJID").alias("donorId"),
                col("AGE").alias("age"),
                col("SEX").alias("sex")
            )
            .withColumn(
                "sex",
                when(col("sex") == "1", lit("M"))
                .when(col("sex") == "2", lit("F"))
                .otherwise(lit("U"))
            )
        )
        
        # Create metadata for each sample column
        # Split OrigSample into donorId + sampleId
        from pyspark.sql.functions import lit as spark_lit
        sample_metadata_rows = []
        for sample_id in self.sample_ids:
            parts = sample_id.split("-", 2)
            donor_id = f"{parts[0]}-{parts[1]}"
            sample_metadata_rows.append((sample_id, donor_id))
        
        # Create dataframe from sample IDs
        df_samples = self.spark.createDataFrame(sample_metadata_rows, ["OrigSample", "donorId"])
        
        # Join with metadata
        df_metadata = (
            df_samples
            .join(meta_sample, on="OrigSample", how="left")
            .join(meta_subject, on="donorId", how="left")
            .withColumn(
                "tissueBiosampleId",
                when(col("TissueOntologyID").isNotNull(), 
                     concat_ws("_", split(col("TissueOntologyID"), ":")))
                .otherwise(lit(None))
            )
            .withColumn("unit", lit("TPM"))
            .withColumn("datasourceId", lit("gtex"))
            .withColumn("datatypeId", lit("bulk rna-seq"))
            .select("OrigSample", "donorId", "tissueBiosampleFromSource", "tissueBiosampleId",
                   "age", "sex", "unit", "datasourceId", "datatypeId")
        )
        
        # Convert to pandas and save as TSV
        print("Converting matrix to pandas...")
        matrix_pd = self.df_matrix.toPandas()
        metadata_pd = df_metadata.toPandas()
        
        # Sort by targetId and sample ID for consistent output
        matrix_pd = matrix_pd.sort_values("targetId")
        metadata_pd = metadata_pd.sort_values("OrigSample")
        
        # Save matrix as compressed TSV
        matrix_output = f"{output_path}matrix/gtex_expression_matrix.tsv.gz"
        metadata_output = f"{output_path}matrix/gtex_sample_metadata.tsv.gz"
        
        # Create directory if it doesn't exist (for local mode)
        import os
        dir_path = f"{output_path}matrix/"
        if dir_path.startswith("file://"):
            dir_path = dir_path.replace("file://", "")
            matrix_output = matrix_output.replace("file://", "")
            metadata_output = metadata_output.replace("file://", "")
        
        os.makedirs(dir_path, exist_ok=True)
        
        print(f"Writing matrix to {matrix_output}...")
        matrix_pd.to_csv(matrix_output, sep="\t", index=False, compression="gzip")
        
        print(f"Writing metadata to {metadata_output}...")
        metadata_pd.to_csv(metadata_output, sep="\t", index=False, compression="gzip")
        
        print(f"Matrix data written to {matrix_output}")
        print(f"Metadata written to {metadata_output}")
        print(f"Matrix shape: {matrix_pd.shape[0]} genes x {matrix_pd.shape[1]-1} samples")

    def main(self):
        self.spark = SparkSession.builder \
                .appName("GTExUnaggregatedExpression") \
                .config("spark.jars", "https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar") \
                .config("spark.executor.memory", "70g") \
                .config("spark.driver.memory", "50g") \
                .config("spark.memory.offHeap.enabled",True) \
                .config("spark.memory.offHeap.size","16g") \
                .getOrCreate()
        print("Reading GTEx data...")
        self.read_gtex_data()
        print("Packing data for output...")
        self.pack_data_for_output(local=self.local, json=self.json, matrix=self.matrix)
        self.spark.stop()

    def __init__(
        self, gtex_source_data_path: str,
        output_directory_path: str,
        sample_metadata_path: str,
        subject_metadata_path: str, 
        json: bool = False,
        local: bool = True,
        matrix: bool = False
    ):
        self.gtex_source_data_path = gtex_source_data_path
        self.output_directory_path = output_directory_path
        self.sample_metadata_path = sample_metadata_path
        self.subject_metadata_path = subject_metadata_path
        self.json = json
        self.local = local
        self.matrix = matrix
        self.spark = None  # Will be initialized in main()

parser = argparse.ArgumentParser(description="Generate unaggregated baseline expression data from GTEx V10.")
parser.add_argument(
    "--local", action="store_true", default=False, 
    help="Run in local mode" 
)
parser.add_argument(
    "--gtex-source-data-path", required=True, type=str, 
    help="GTEx downloaded sample-by-gene TPM counts, in gzipped GCT format."
)
parser.add_argument(
    "--output-directory-path", required=True, type=str, default='.',
    help="A path to output the output file. One object per line.",
)
parser.add_argument(
    "--sample-metadata-path", required=True, type=str, help="A path to the sample metadata file.",
)
parser.add_argument(
    "--subject-metadata-path", required=True, type=str, help="A path to the subject metadata file.",
)
parser.add_argument(
    "--json", action='store_true', default=False,
     help="Save output as JSON instead of parquet.",
)
parser.add_argument(
    "--matrix", action='store_true', default=False,
     help="Save output as matrix form (genes x samples) with separate metadata TSV files, both gzipped.",
)


if __name__ == "__main__":
    args = parser.parse_args()
    BaselineExpression(
        args.gtex_source_data_path,
        args.output_directory_path,
        args.sample_metadata_path,
        args.subject_metadata_path,
        args.json,
        args.local,
        args.matrix
    ).main()
