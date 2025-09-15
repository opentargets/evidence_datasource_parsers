from pyspark.sql import SparkSession
from pyspark.sql import functions as f
from functools import reduce
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

    def override_biosample_id(self, biosample_mapping_reference_path: str,
                                key_col: str,
                                source_col: str,
                                id_col: str):
        """
        Assign a new biosample Id to each row in self.df by exact matching
        `source_col` to entries in a wide TSV mapping file.

        Parameters
        ----------
        biosample_mapping_reference_path : str
            Path to the wide mapping file (TSV with header).
        key_col : str
            The column name in the mapping file that holds the canonical BiosampleId (e.g., "BiosampleId").
        source_col : str
            The column name in self.df with the tissue string to match (e.g., "tissueBiosampleFromSource").
        id_col : str
            The column name in self.df to overwrite when a match is found (e.g., "tissueBiosampleId").

        Returns
        -------
        pyspark.sql.DataFrame
            Updated DataFrame (also stored back into self.df).
        """

        spark = self.spark
        df = self.df  # your expression dataframe

        # ---- checks ----
        if source_col not in df.columns or id_col not in df.columns:
            raise ValueError(f"`self.df` must contain '{source_col}' and '{id_col}'")

        # ---- load mapping (TSV) ----
        mapping_df = (
            spark.read
            .option("header", True)
            .option("inferSchema", True)
            .option("sep", "\t")
            .csv(biosample_mapping_reference_path)
        ).drop('Count')  # drop any extraneous columns
        if key_col not in mapping_df.columns:
            raise ValueError(f"Mapping file must contain '{key_col}'")

        other_cols = [c for c in mapping_df.columns if c != key_col]
        if not other_cols:
            raise ValueError("Mapping file must have at least one non-key column with tokens to match")

        # ---- helpers (UDF-free) ----
        def norm(col):
            # Exact-match key: lowercased, non-alphanumerics removed.
            return f.regexp_replace(f.lower(f.col(col)), r"[^a-z0-9]", "")

        def pretty(col):
            """
            Human-friendly display:
            - '_'/'-' → spaces
            - split camelCase and 'TReg' → 'T Reg'
            - keep acronyms/digit mixes (CD4, HLA, NK) and single-letter caps as-is
            - lowercase other words
            - remove , and .
            """
            s = f.coalesce(f.col(col).cast("string"), f.lit(""))

            # unify separators
            s = f.regexp_replace(s, r"[_-]+", " ")
            # 'TReg' -> 'T Reg' (without breaking 'CD4')
            s = f.regexp_replace(s, r"([A-Z])([A-Z][a-z])", r"\1 \2")
            # split camelCase: adrenalGland -> adrenal Gland
            s = f.regexp_replace(s, r"(?<=[a-z])(?=[A-Z])", " ")
            # number/letter boundaries (won't split inside 'CD4')
            s = f.regexp_replace(s, r"(?<=[a-z])(?=\d)", " ")
            s = f.regexp_replace(s, r"(?<=\d)(?=[a-z])", " ")
            # remove , and .
            s = f.regexp_replace(s, r"[.,]", "")
            # normalize spaces
            s = f.regexp_replace(s, r"\s+", " ")
            s = f.trim(s)

            words = f.split(s, " ")
            words2 = f.transform(
                words,
                lambda w: f.when(
                    w.rlike(r'^[A-Z0-9]{2,}$') | w.rlike(r'^[A-Z]$'),  # keep CD4/HLA/NK/T/B
                    w
                ).otherwise(f.lower(w))
            )
            return f.array_join(words2, " ")

        def clean_tokens(colname, delim=";"):
            s = f.coalesce(f.col(colname).cast("string"), f.lit(""))
            arr = f.split(s, delim)
            arr = f.transform(arr, lambda x: f.trim(x))
            return f.filter(arr, lambda x: x.isNotNull() & (x != "") & (f.lower(x) != f.lit("null")))

        # ---- build long mapping directly (skip any list middleman) ----
        long_parts = [
            mapping_df.select(f.col(key_col), f.explode_outer(clean_tokens(c)).alias("value"))
            for c in other_cols
        ]
        mapping_long = reduce(lambda a, b: a.unionByName(b), long_parts) \
            .filter(f.col("value").isNotNull() & (f.col("value") != ""))

        # Enrich with norm & pretty; choose ONE row per norm deterministically (min key).
        enriched = (
            mapping_long
            .withColumn("norm", norm("value"))
            .withColumn("pretty", pretty("value"))
            .filter(f.col("norm") != "")
        )
        min_key_by_norm = enriched.groupBy("norm").agg(f.min(f.col(key_col)).alias(key_col))
        mapping_norm = (
            enriched.join(min_key_by_norm, on=["norm", key_col], how="inner")
            .select("norm", key_col, "pretty")
            .dropDuplicates(["norm"])  # one mapping per normalized term
        )

        # ---- normalize df[source_col], join, and replace ----
        df_norm = df.withColumn("norm", norm(source_col)).drop(id_col)
        joined = df_norm.join(mapping_norm, on="norm", how="left")

        df_updated = (
            joined
            .withColumn(id_col, f.col(key_col))
            .withColumn(source_col, f.coalesce(f.col("pretty"), f.col(source_col)))
            .drop("norm", key_col, "pretty")
        )

        # store and return
        self.df = df_updated

    def drop_null_biosample_ids(self):
        """
        Drop rows where any of the specified celltypeBiosampleId AND tissueBiosampleId are null
        """
        # First check if the columns exist and add them if they don't, populating with nulls
        if "tissueBiosampleId" not in self.df.columns:
            self.df = self.df.withColumn("tissueBiosampleId", f.lit(None).cast("string"))
            self.df = self.df.withColumn("tissueBiosampleFromSource", f.lit(None).cast("string"))
        if "celltypeBiosampleId" not in self.df.columns:
            self.df = self.df.withColumn("celltypeBiosampleId", f.lit(None).cast("string"))
            self.df = self.df.withColumn("celltypeBiosampleFromSource", f.lit(None).cast("string"))

        # Find the rows where both tissue and celltype biosample ID are null
        nulls = self.df.filter(
            (self.df["tissueBiosampleId"].isNull()) & (self.df["celltypeBiosampleId"].isNull())
        )
        print("The following biosample from source have null biosample IDs:")
        print(nulls.select("tissueBiosampleFromSource", "celltypeBiosampleFromSource").distinct().show(100,truncate=False))
        # Then drop those rows from the dataframe
        self.df = self.df.filter(~((self.df["tissueBiosampleId"].isNull()) & (self.df["celltypeBiosampleId"].isNull())))

    def within_donor_mean(self,
                        tissue_col: str = "tissueBiosampleFromSource",
                        celltype_col: str = "celltypeBiosampleFromSource",
                        expr_col: str = "expression",
                        sep: str = ", "):
        """
        Group by ALL columns except `expr_col`, `tissue_col`, and `celltype_col`.

        After grouping:
        - `expr_col` becomes the mean of expressions in the group
        - `tissue_col` becomes a comma-joined concatenation of all non-null/non-empty
            entries in the group, sorted lexicographically (duplicates preserved)
        - `celltype_col` same as above

        Result is stored back to self.df
        """

        df = self.df

        # Columns to group by
        group_cols = [c for c in df.columns if c not in {expr_col, tissue_col, celltype_col}]

        if group_cols:
            agg_df = (
                df.groupBy(group_cols)
                .agg(
                    f.avg(f.col(expr_col)).alias(expr_col),
                    f.collect_set(f.col(tissue_col)).alias("_tissues"),
                    f.collect_set(f.col(celltype_col)).alias("_celltypes"),
                )
            )
        else:
            # If nothing to group by, aggregate globally
            agg_df = df.agg(
                f.avg(f.col(expr_col)).alias(expr_col),
                f.collect_set(f.col(tissue_col)).alias("_tissues"),
                f.collect_set(f.col(celltype_col)).alias("_celltypes"),
            )

        # Helper: trim -> drop null/empty -> sort -> join
        def join_sorted(arr_col, out_name):
            arr_trim = f.transform(f.col(arr_col), lambda x: f.when(x.isNull(), None).otherwise(f.trim(x)))
            arr_keep = f.filter(arr_trim, lambda x: x.isNotNull() & (x != ""))
            arr_sorted = f.array_sort(arr_keep)  # simple lexicographic sort; duplicates preserved
            return f.when(f.size(arr_sorted) > 0, f.array_join(arr_sorted, sep)).otherwise(f.lit(None)).alias(out_name)

        agg_df = (
            agg_df
            .withColumn(tissue_col, join_sorted("_tissues", tissue_col))
            .withColumn(celltype_col, join_sorted("_celltypes", celltype_col))
            .drop("_tissues", "_celltypes")
        )

        self.df = agg_df

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
        self.df = quartile_df

    def calculate_expression_distribution(self, local=False):
        """
        This function groups the dataframe by datasourceId, targetId and datatypeId
        then calculates the distribution of expression values for each gene e.g.
        if a gene is expressed (> 0) in 7 out of a possible 10 samples, its distribution would be 0.7.
        It works on the median expression value.
        """
        # Group by the relevant columns and count the number of non-zero expressions
        exp_distribution_df = self.df.groupBy(
            "targetId",
            "datasourceId",
            "datatypeId",
            "unit"
        ).agg(
            (f.sum(f.when(f.col("median") > 0, 1).otherwise(0)) / f.count("*")).alias("distribution")
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
    "--override-before", action='store_true', default=False,
    help="Flag indicating to override before aggregating expression data rather than after.",
)
parser.add_argument(
    "--override-tissue-reference", type=str, default=None,
    help="A path to a TSV containing which annotations should be overridden ",
)
parser.add_argument(
    "--override-celltype-reference", type=str, default=None,
    help="A path to a TSV containing which annotations should be overridden ",
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
    override_before = args.override_before
    override_tissue_reference = args.override_tissue_reference
    override_celltype_reference = args.override_celltype_reference

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

    if override_before:
        if override_tissue_reference is not None:
            print("Overriding tissue biosample IDs before aggregation...")
            eq.override_biosample_id(
                biosample_mapping_reference_path=override_tissue_reference,
                key_col="BiosampleId",
                source_col="tissueBiosampleFromSource",
                id_col="tissueBiosampleId"
            )
        if override_celltype_reference is not None:
            print("Overriding celltype biosample IDs before aggregation...")
            eq.override_biosample_id(
                biosample_mapping_reference_path=override_celltype_reference,
                key_col="BiosampleId",
                source_col="celltypeBiosampleFromSource",
                id_col="celltypeBiosampleId"
            )
        print("Dropping rows where both tissue AND celltype biosample ID are null...")
        eq.drop_null_biosample_ids()

        print("Calculating within-donor mean expression...")
        eq.within_donor_mean()

    print("Calculating quartiles...")
    quartile_df = eq.calculate_quartiles(local=local)

    if not override_before:
        if override_tissue_reference is not None:
            print("Overriding tissue biosample IDs after aggregation...")
            eq.override_biosample_id(
                biosample_mapping_reference_path=override_tissue_reference,
                key_col="BiosampleId",
                source_col="tissueBiosampleFromSource",
                id_col="tissueBiosampleId"
            )
        if override_celltype_reference is not None:
            print("Overriding celltype biosample IDs after aggregation...")
            eq.override_biosample_id(
                biosample_mapping_reference_path=override_celltype_reference,
                key_col="BiosampleId",
                source_col="celltypeBiosampleFromSource",
                id_col="celltypeBiosampleId"
            )
        print("Dropping rows where both tissue AND celltype biosample ID are null...")
        eq.drop_null_biosample_ids()
    
    print("Calculating expression distribution...")
    distribution_df = eq.calculate_expression_distribution()

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