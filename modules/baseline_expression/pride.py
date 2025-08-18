"""Ingests PRIDE proteomic data and generates the unaggregated baseline expression data.

Requires Apache Spark (pyspark) for distributed processing.

Expected input formats:
- PRIDE source data: TSV file containing sample-by-protein PPB counts.
- Sample metadata: JSON file with experimental designs.
- Tissue to ontology mapping: TSV file with columns PROPERTY VALUE, ONTOLOGY TERM(S)
- Target index: Parquet file with target IDs and associated protein IDs.
"""

import argparse
import functools
import gzip
import json
import os

from pyspark.sql.functions import split, col, array, struct, lit, explode, concat_ws, when, regexp_replace, element_at, concat
from pyspark.sql import SparkSession

class BaselineExpression:
    """Collection of steps to generate unaggregated baseline expression data."""

    def read_pride_data(
        self,
        pride_code: str
    ):
        pride_source_data_path = f"{self.pride_source_data_dir}/{pride_code}/{pride_code}_OpenTargets_Quant_ppb.txt"
        pride_sdrf_path = f"{self.pride_source_data_dir}/{pride_code}/{pride_code}_OpenTargets_SDRF.json"
        # Read the PRIDE matrix file
        pride_matrix = (
            self.spark.read
            .option("sep", "\t")
            .option("header", "true")
            .csv(pride_source_data_path)
        )
        
        # Extract the uniprot IDs and bind in the target index for Ensembl IDs
        target_index = self.spark.read.parquet(self.target_index_path)
        target_mapping = target_index.select(
            col("id"),
            explode(col("proteinIds")).alias("proteinId")
        ).select(
            col("id"),
            col("proteinId.id").alias("proteinId"),
        )
        
        # Extract the uniprot IDs and bind in the target index for Ensembl IDs
        orig_cols = pride_matrix.columns

        pride_matrix = (
            pride_matrix
            .withColumn(
                "proteinId",
                element_at(split(col("Protein IDs"), r"\|"), 2)
            )
            .withColumn("proteinId", regexp_replace("proteinId", r"-\d+$", ""))  # optional, removes isoform suffix
            .join(target_mapping, on="proteinId", how="left")
            .select("id", "proteinId", "Gene Symbol", "Protein IDs", *orig_cols)
            .drop("ENSG")
        )

        # Wide→long: explode an array of structs (one per sample column)
        data_cols = [c for c in pride_matrix.columns if c not in ["id","proteinId","Gene Symbol","Protein IDs"]]

        # build an array<struct<Sample:string,PPB:double>>  
        samp_structs = array(*[
            struct(lit(c).alias("Sample"), col(c).cast("double").alias("PPB"))
            for c in data_cols
        ]).alias("samp_ppb")

        pride_long = (
            pride_matrix
            .select("id", "proteinId","Gene Symbol", "Protein IDs", explode(samp_structs).alias("x"))
            .select(
                col("id"),
                col("proteinId"),
                col("Gene Symbol"),
                col("Protein IDs"),
                col("x.Sample").alias("assayId"),
                col("x.PPB")
            )
        )

        # Read the PRIDE SDRF file
        pride_sdrf = (
            self.spark.read
            .option("multiline", True)   # important for JSON
            .json(pride_sdrf_path)
            .withColumn("ex", explode("experimentalDesigns"))
            .select("experimentId", col("ex.*"))
            # In the sex column, replace not available' with NULL
            .withColumn("sex", when(col("sex") == "not available", None).otherwise(col("sex"))) 
            .withColumn("age", when(col("age") == "not available", None).otherwise(col("age")))
        )

        # Read the PRIDE ontology mapping
        pride_ontology_mapping = (
            self.spark.read
            .option("header", "true")
            .option("sep", "\t")
            .csv(self.tissue_ontology_mapping_path)
            .withColumnRenamed("PROPERTY VALUE", "tissue")
            .withColumnRenamed("ONTOLOGY TERM(S)", "tissueOntologyTerm")
            .select("tissue", "tissueOntologyTerm") # Only select the relevant terms
        )

        # Join the long DataFrame with the SDRF metadata and ontology mapping
        pride_long = (
            pride_long
            .join(pride_sdrf, on="assayId", how="left")
            .withColumn("donorId", concat(col("experimentId"), lit("-"), col("individual"))) # Add a column for donor
            .join(pride_ontology_mapping, on="tissue", how="left")
            .select(
                "id", "proteinId", "donorId", "PPB", "tissueOntologyTerm",
                "age", "sex", "tissue"
            )
        )

        # Drop rows where tissueOntologyTerm is null
        pride_long = pride_long.filter(col("tissueOntologyTerm").isNotNull())

        pride_long = (pride_long
            .withColumnRenamed("id", "targetId")
            .withColumnRenamed("PPB", "expression")
            .withColumnRenamed("tissue", "tissueBiosampleFromSource")
            .withColumnRenamed("tissueOntologyTerm", "tissueBiosampleId")
            .withColumnRenamed("proteinId", "targetFromSourceId")
            .withColumn("unit", lit("PPB (iBAQ)"))  # Add unit column with constant value "PPB (iBAQ)"
            .withColumn("datatypeId", lit("mass-spectrometry proteomics"))  # Add datatypeId column with constant value "mass-spectrometry proteomics"
            .withColumn("datasourceId", lit("PRIDE"))  # Add datasourceId column with constant value "PRIDE"
        )
        # Bind the pride long dataframe to the existing self.df DataFrame
        if hasattr(self, 'df'):
            self.df = self.df.unionByName(pride_long, allowMissingColumns=True)
        else:
            self.df = pride_long


    def pack_data_for_output(self, local: bool = False, json: bool = False):
        """Use spark to write the DataFrame to parquet format."""
        if local:
            output_path = f"file://{self.output_directory_path}"
        else:
            output_path = f"{self.output_directory_path}"
        if json:
            output_path = f"{output_path}/json/pride_baseline_expression"
            # If JSON output is requested, convert DataFrame to JSON format
            self.df.write.mode("overwrite").json(output_path)
            print(f"Data written to {output_path} in JSON format")
        else:
            output_path = f"{output_path}/parquet/pride_baseline_expression"
            # If parquet output is requested, convert DataFrame to parquet format
            self.df.write.mode("overwrite").parquet(output_path)
            print(f"Data written to {output_path} in Parquet format")


    def main(self):
        self.spark = SparkSession.builder \
                .appName("PRIDEUnaggregatedExpression") \
                .config("spark.jars", "https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar") \
                .config("spark.executor.memory", "70g") \
                .config("spark.driver.memory", "50g") \
                .config("spark.memory.offHeap.enabled",True) \
                .config("spark.memory.offHeap.size","16g") \
                .getOrCreate()
        print("Reading PRIDE data...")
        for pride_code in self.pride_codes:
            print(f"Processing PRIDE code: {pride_code}")
            self.read_pride_data(pride_code)
        print("Packing data for output...")
        self.pack_data_for_output(local=self.local, json=self.json)
        self.spark.stop()

    def __init__(
        self,
        pride_source_data_dir: str,
        pride_codes: str,
        output_directory_path: str,
        tissue_ontology_mapping_path: str,
        target_index_path: str,
        json: bool = False,
        local: bool = True
    ):
        self.pride_source_data_dir = pride_source_data_dir
        self.pride_codes = pride_codes.split(",")  # Split multiple codes if provided
        self.output_directory_path = output_directory_path
        self.tissue_ontology_mapping_path = tissue_ontology_mapping_path
        self.target_index_path = target_index_path
        self.json = json
        self.local = local



parser = argparse.ArgumentParser(description="Generate unaggregated baseline expression data from PRIDE.")
parser.add_argument(
    "--local", action="store_true", default=False, 
    help="Run in local mode" 
)
parser.add_argument(
    "--pride-source-data-dir", required=True, type=str,
    help="Directory containing PRIDE data and metadata files in subdirectories."
)
parser.add_argument(
    "--pride-codes", required=True, type=str,
    help="PRIDE codes for the datasets e.g. 'PXD000000' or 'PXD000000,PXD000001'.",
)
parser.add_argument(
    "--output-directory-path", required=True, type=str, default='.',
    help="A path to output the output file. One object per line.",
)
parser.add_argument(
    "--tissue-ontology-mapping", required=True, type=str, help="A path to the tissue ontology mapping file.",
)
parser.add_argument(
    "--json", action='store_true', default=False,
     help="Save output as JSON instead of parquet.",
)
parser.add_argument(
    "--target-index", required=True, type=str,
    help="A path to the target index parquet file.",
)


if __name__ == "__main__":
    args = parser.parse_args()
    BaselineExpression(
        args.pride_source_data_dir,
        args.pride_codes,
        args.output_directory_path,
        args.tissue_ontology_mapping,
        args.target_index,
        args.json,
        args.local
    ).main()


# Run the script with:
# python pride.py --pride-codes PXD000000,PXD000001 --pride-source-data-dir /path/to/pride/data  --output-directory-path /path/to/output --tissue-ontology-mapping /path/to/tissue/ontology/mapping.tsv
# python /home/alegbe/repos/evidence_datasource_parsers/modules/baseline_expression/pride.py --pride-codes PXD020192,PXD000865,PXD000561 --pride-source-data-dir /home/alegbe/data/bulkProteomics/ --output-directory-path /home/alegbe/results/unaggregated/pride/ --tissue-ontology-mapping /home/alegbe/data/bulkProteomics/pride_ontology_mapping_sorted.txt --target-index gs://open-targets-ppp-data-releases/25.06/output/target/
