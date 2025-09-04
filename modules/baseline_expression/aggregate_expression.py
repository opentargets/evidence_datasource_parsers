from pyspark.sql import SparkSession
from pyspark.sql import functions as f
import os
import argparse
import pandas as pd

class AggregateExpression:
    """
    This class is used to take the average of the Pseudobulk expression data so we have a single value per gene per annotation. The output of this class is a dataframe with the average expression of each gene (rows) per annotation (columns)
    """

    def load_data(self, directory, local=False):
        """
        This function loads the data from the directory and returns a list of DataFrames.
        """
        # Use spark to read the parquet files in the directory
        if local:
            directory = f"file://{directory}"
        df = self.spark.read.parquet(directory)

        self.df = df

    def within_donor_mean(self, annotation_groupings):
        """
        This function calculates the mean expression of each gene within each grouping
        of datasourceId, tissueBiosampleFromSource, celltypeBiosampleFromSource and donorId.
        """

    def calculate_quartiles(self, local=False):
        """
        This function calculates the expression quartiles of each gene across all donors.
        """
       # Define the grouping cols
        groupby_cols = [
            "targetId",
            "datasourceId",
            "datatypeId",
            "unit"
        ]
        if 'tissueBiosampleFromSource' in self.df.columns:
            groupby_cols.append("tissueBiosampleFromSource")
            groupby_cols.append("tissueBiosampleId")
        if 'celltypeBiosampleFromSource' in self.df.columns:
            groupby_cols.append("celltypeBiosampleFromSource")
            groupby_cols.append("celltypeBiosampleId")
        
        # Partition by grouping keys
        quartile_df = self.df.repartition(
            *groupby_cols
        )

        # Define the quantile probabilities
        quartile_probs = [0,0.25, 0.50, 0.75, 1]

        # Group and compute approximate quantiles, default params
        quartile_df = quartile_df.groupBy(
            *groupby_cols
        ).agg(
            f.percentile_approx("expression", quartile_probs).alias("q_vals"))


        quartile_df = quartile_df.select(
            *groupby_cols,
            f.col("q_vals")[0].alias("min"),
            f.col("q_vals")[1].alias("q1"),
            f.col("q_vals")[2].alias("median"),
            f.col("q_vals")[3].alias("q3"),
            f.col("q_vals")[4].alias("max")
        )
        return quartile_df

    def calculate_expression_distribution(self, local=False):
        """
        This function groups the dataframe by datasourceId, targetId and datatypeId
        then calculates the distribution of expression values for each gene e.g.
        if a gene is expressed (> 0) in 7 out of a possible 10 samples, its distribution would be 0.7.
        It works on the maximum expression value.
        """
        # Group by the relevant columns and count the number of non-zero expressions
        exp_distribution_df = self.df.groupBy(
            "targetId",
            "datasourceId",
            "datatypeId",
            "unit"
        ).agg(
            (f.sum(f.when(f.col("max") > 0, 1).otherwise(0)) / f.count("*")).alias("distribution")
        )
        return exp_distribution_df

    def write_data(self, output_directory, json = False):
        """
        This function writes the DataFrame to parquet format.
        """
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        if json:
            self.df.write.mode("overwrite") \
                .json(output_directory)
        else:
            # If not JSON, write as parquet
            self.df.write.mode("overwrite") \
                .parquet(output_directory)
            

    def __init__(self, local=False):
        """
        This function initializes the class.
        """
        if local:
            self.spark = SparkSession.builder.master("local").appName("spark_etl").config("spark.hadoop.fs.defaultFS", "file:///").getOrCreate()
        else:
            self.spark = SparkSession.builder.appName("ProcessExpression").getOrCreate()
        
        # Disable whole-stage code generation
        self.spark.conf.set("spark.sql.codegen.wholeStage", "false")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--local", action="store_true", default=False, 
    help="Run in local mode" 
)
parser.add_argument(
    "--directory", required=True, type=str, 
    help="Directory containing expression data"
)
parser.add_argument(
    "--json", action='store_true', default=False,
    help="Save output as JSON instead of parquet.",
)
parser.add_argument(
    "--combined-annots-reference", action='store_true', default=False,
    help="A dataframe containing which annotations should be combined ",
)
parser.add_argument(
    "--output", required=True, type=str, 
    help="Output directory"
)

if __name__ == "__main__":
    args = parser.parse_args()
    directory = args.directory
    local = args.local
    json = args.json
    output_dir = args.output

    eq = AggregateExpression(local=local)
    eq.spark = SparkSession.builder \
        .appName("AggregateExpression") \
        .config("spark.jars", "https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar") \
        .config("spark.executor.memory", "70g") \
        .config("spark.driver.memory", "50g") \
        .config("spark.memory.offHeap.enabled", True) \
        .config("spark.memory.offHeap.size", "16g") \
        .getOrCreate()
    eq.load_data(directory, local=local)
    print("Calculating expression distribution...")
    distribution_df = eq.calculate_expression_distribution()
    print("Calculating quartiles...")
    quartile_df = eq.calculate_quartiles(local=local)
    print("Packing data for output...")
    file_format = "json" if json else "parquet"
    # Pull out the output directory structure from the input directory
    parts = [seg for seg in os.path.normpath(directory).split(os.sep) if seg]
    datasource = parts[-3]
    dataset_name   = parts[-1]
    output_dir = f"{output_dir}/{datasource}/{file_format}/{dataset_name}"
    # Write the data to the output directorys
    eq.write_data(output_dir, json=json)
    print(f"Data written to {output_dir}")
    eq.spark.stop()