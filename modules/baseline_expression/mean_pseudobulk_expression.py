from pyspark.sql import SparkSession
from pyspark.sql.functions import col, expr
from functools import reduce
import os
import argparse

class AveragePseudobulkExpression:
    """
    This class is used to take the average of the Psuedobulk expression data so we have a single value per gene per annotation. The output of this class is a dataframe with the average expression of each gene (rows) per annotation (columns)
    """
    def list_files(self, directory):
        """
        This function lists all the files in a directory.
        """
        return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".tsv")]

    def calculate_average(self, local=False):
        """
        This function calculates the average expression of each gene across all donors.
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
            
            # Compute the average across replicate columns.
            # This builds an expression like: (`col1` + `col2` + ...)/N.
            avg_expr = "(" + " + ".join([f"`{c}`" for c in donor_cols]) + f")/{len(donor_cols)}"
            
            # Select gene id and the computed average (renaming the column to the cell type).
            df_avg = df.select("ensg", expr(avg_expr).alias(annotation))

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
            self.spark = SparkSession.builder.appName("AveragePseudobulkExpression").getOrCreate()

parser = argparse.ArgumentParser()
parser.add_argument(
    "--local", action="store_true", help="Run in local mode", default=False
)
parser.add_argument(
    "--directory", required=True, type=str, help="Directory containing pseudobulk expression data"
)
parser.add_argument(
    "--output", required=True, type=str, help="Output directory"
)

if __name__ == "__main__":
    args = parser.parse_args()
    directory = args.directory
    local = args.local
    output_dir = args.output

    ape = AveragePseudobulkExpression()
    file_list = ape.list_files(directory)
    dfs = []
    ape.calculate_average()
    df = reduce(lambda x, y: x.join(y, on="ID", how="outer"), dfs)
    outpath = os.path.join(output_dir, "average_pseudobulk_expression.tsv")
    if local:
        df.write.option("header", "true").option("sep", "\t").csv(f"file://{outpath}")
    else:
        df.write.option("header", "true").option("sep", "\t").csv("average_pseudobulk_expression.tsv")
    ape.spark.stop()