from pyspark.sql import SparkSession
from pyspark.sql.functions import col, expr
from functools import reduce
import os
import argparse
import pandas as pd

class ExpressionQuartiles:
    """
    This class is used to take the average of the Pseudobulk expression data so we have a single value per gene per annotation. The output of this class is a dataframe with the average expression of each gene (rows) per annotation (columns)
    """
    def list_files(self, directory):
        """
        This function lists all the files in a directory.
        """
        return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".tsv")]

    def calculate_quartiles(self, local=False):
        """
        This function calculates the expression quartiles of each gene across all donors.
        """
        for i,file in enumerate(file_list):
            print(f"Processing file {i+1}/{len(file_list)}: {file}")
            # try:
            # Read the file as a DataFrame (tab-separated with header)
            if local:
                df = self.spark.read.option("header", "true").option("sep", "\t").csv(f"file://{file}")
            else:
                df = self.spark.read.option("header", "true").option("sep", "\t").csv(file)
            
            # Identify donor columns (all columns except "ensg")
            donor_cols = [d for d in df.columns if d != "ensg"]
            
            # Cast donor columns to Double.
            for donor in donor_cols:
                df = df.withColumn(donor, col(donor).cast("double"))
            
            # Extract the annotation from the file name.
            annotation = os.path.basename(file).split(".")[0]
            
            # Use the ApproxQuantile function to compute the quartiles.

            quartiles = df.approxQuantile(donor_cols, [0, 0.25, 0.5, 0.75, 1], 0.01)
            # Create a new DataFrame with the quartiles.
            df_quartiles = self.spark.createDataFrame(
                [(annotation, *quartiles)],
                schema=["annotation", "Min", "Q1", "Q2", "Q3", "Max"]
            )
            

            # Rename gene id from "ensg" to "ID".
            df_avg = df_avg.withColumnRenamed("ensg", "ID")
            
            dfs.append(df_avg)
            # except Exception as e:
            #     print(f"Error processing {file}: {e}")

    def __init__(self):
        """
        This function initializes the class.
        """
        if local:
            self.spark = spark = SparkSession.builder.master("local").appName("spark_etl").config("spark.hadoop.fs.defaultFS", "file:///").getOrCreate()
        else:
            self.spark = SparkSession.builder.appName("ProcessExpression").getOrCreate()
        
        # Disable whole-stage code generation
        self.spark.conf.set("spark.sql.codegen.wholeStage", "false")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--local", action="store_true", help="Run in local mode", default=False
)
parser.add_argument(
    "--directory", required=True, type=str, help="Directory containing expression data"
)
parser.add_argument(
    "--schema", required=True, type=str, help="Path to the expression_aggregated.json schema"
)
parser.add_argument(
    "--output", required=True, type=str, help="Output directory"
)

if __name__ == "__main__":
    args = parser.parse_args()
    directory = args.directory
    local = args.local
    output_dir = args.output

    eq = ExpressionQuartiles()
    file_list = eq.list_files(directory)
    dfs = []
    eq.calculate_quartiles()
    
    df_list = [df_avg.toPandas() for df_avg in dfs]
    print('Merging dataframes')
    df_merged = reduce(lambda left, right: pd.merge(left, right, on="ID", how="outer"), df_list)
    # outpath = os.path.join(output_dir, "average_pseudobulk_expression.tsv")
    print(f"Writing to {outpath}")
    if local:
        df_merged.to_csv(outpath, sep="\t", index=False)
    else:
        df_merged.to_csv(f"file://{outpath}", sep="\t", index=False)
    eq.spark.stop()