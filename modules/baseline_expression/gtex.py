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

from pyspark.sql.functions import split, col, array, struct, lit, explode, concat_ws, when


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
            .withColumn("Name", split(col("Name"), r"\.").getItem(0))  # strip .2 etc.
        )

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
                col("SAMPID").alias("sampleId"),
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
        self.df = (
            df_long
            .join(meta_sample,  on="sampleId", how="left")
            .join(meta_subject, on="donorId",    how="left")
            .select(
                "Name", "donorId", "sampleId", "TPM",
                "Tissue", "TissueOntologyID", "Age", "Sex"
            )
        )

    def pack_data_for_output(self, local: bool = False):
        """Use spark to write the DataFrame to JSONL format."""
        import os

        import os
        output_path = os.path.abspath(os.path.join(self.output_directory_path, "gtex_baseline_expression"))
        self.df.write.mode("overwrite").json(output_path)
        print(f"Data written to {output_path}")
        # If local mode, write to a local file
        self.df.write.mode("overwrite").json(f"file://{self.output_directory_path}/gtex_baseline_expression")
        print(f"Data written to {self.output_directory_path}/gtex_baseline_expression")

    def main(self):
        from pyspark.sql import SparkSession
        self.spark = SparkSession.builder.appName("GTExBaselineExpression").getOrCreate()
        print("Reading GTEx data...")
        self.read_gtex_data()
        print("Packing data for output...")
        self.pack_data_for_output(local=True)
        self.spark.stop()

    def __init__(
        self, gtex_source_data_path: str,
        output_directory_path: str,
        sample_metadata_path: str,
        subject_metadata_path: str, 
    ):
        self.gtex_source_data_path = gtex_source_data_path
        self.output_directory_path = output_directory_path
        self.sample_metadata_path = sample_metadata_path
        self.subject_metadata_path = subject_metadata_path
        self.spark = None  # Will be initialized in main()

parser = argparse.ArgumentParser(description="Generate unaggregated baseline expression data from GTEx V10.")
parser.add_argument(
    "--gtex-source-data-path", required=True, type=str, 
    help="GTEx downloaded sample-by-gene TPM counts, in gzipped GCT format."
)
parser.add_argument(
    "--output-directory-path", required=True, type=str, default='.',
    help="A path to output the output file. GZIP-compressed JSON, one object per line.",
)
parser.add_argument(
    "--sample-metadata-path", required=True, type=str, help="A path to the sample metadata file.",
)
parser.add_argument(
    "--subject-metadata-path", required=True, type=str, help="A path to the subject metadata file.",
)


if __name__ == "__main__":
    args = parser.parse_args()
    BaselineExpression(
        args.gtex_source_data_path,
        args.output_directory_path,
        args.sample_metadata_path,
        args.subject_metadata_path,
    ).main()
