import argparse
import logging
import sys

from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.types import FloatType
from pyspark.sql.functions import col, collect_set, split, element_at, struct, lit

# A few genes do not have Ensembl IDs in the data file provided
CRISPR_SYMBOL_MAPPING = {
    'CASC5': 'ENSG00000137812',
    'CIRH1A': 'ENSG00000141076',
    'EFTUD1': 'ENSG00000140598',
    'ENSG00000163660': 'ENSG00000163660',
    'KIAA0947': 'ENSG00000164151',
    'KIAA1432': 'ENSG00000107036',
    'NDNL2': 'ENSG00000185115',
    'SRPR': 'ENSG00000182934',
    'ZNF259': 'ENSG00000109917'
}


def main(desc_file, evid_file, cell_file, out_file):
    sparkConf = (
        SparkConf()
        .set('spark.driver.memory', '15g')
        .set('spark.executor.memory', '15g')
        .set('spark.driver.maxResultSize', '0')
        .set('spark.debug.maxToStringFields', '2000')
        .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
    )
    spark = (
        SparkSession.builder
        .config(conf=sparkConf)
        .master('local[*]')
        .getOrCreate()
    )

    # Log parameters:
    logging.info(f'Evidence file: {evid_file}')
    logging.info(f'Description file: {desc_file}')
    logging.info(f'Cell type annotation: {cell_file}')
    logging.info(f'Output file: {out_file}')

    # Read files:
    evidence_df = (
        spark.read.csv(evid_file, sep='\t', header=True)
        .drop('pmid', 'gene_set_name', 'disease_name')
    )
    cell_lines_df = spark.read.csv(cell_file, sep='\t', header=True)
    description_df = spark.read.csv(desc_file, sep='\t', header=True)

    # Logging dataframe stats:
    logging.info(f'Number of evidence: {evidence_df.count()}')
    logging.info(f'Number of descriptions: {description_df.count()}')
    logging.info(f'Number of cell/tissue annotation: {cell_lines_df.count()}')

    # Tissues and cancer types are annotated together in the same column (tissue_or_cancer_type)
    # To disambiguate one from another, the column is combined with the cell lines
    # First on the tissue level:
    tissue_desc = (
        description_df
        .withColumnRenamed('tissue_or_cancer_type', 'tissue')
        .join(cell_lines_df, on='tissue', how='inner')
    )

    # And then on the disease level:
    cell_desc = (
        description_df
        .withColumnRenamed('tissue_or_cancer_type', 'diseaseFromSource')
        .join(cell_lines_df, on='diseaseFromSource', how='inner')
    )

    merged_annotation = (
        # Concatenating the above generated dataframes:
        cell_desc
        .union(tissue_desc)

        # Aggregating by disease and method:
        .groupBy('diseaseFromSource', 'efo_id', 'method')

        # The cell annotation is aggregated in a list of struct:
        .agg(
            collect_set(struct(
                col('name'), col('id'), col('tissue'), col('tissueId')
            )).alias('diseaseCellLines')
        )
        .drop('method')
    )

    # Joining merged annotation with evidence:
    pooled_evidence_df = (
        evidence_df
        .select(
            col('target_id').alias('targetFromSourceId'),
            col('disease_id').alias('efo_id'),
            col('score').alias('resourceScore').cast(FloatType()),
        )

        # Some of the target identifier are not Ensembl Gene id - replace them:
        .replace(to_replace=CRISPR_SYMBOL_MAPPING, subset=['target_id'])

        # Merging with descriptions:
        .join(merged_annotation, on='efo_id', how='outer')

        # From EFO uri, generate EFO id:
        .withColumn(
            'diseaseFromSourceMappedId',
            element_at(split(col('efo_id'), '/'), -1).alias('diseaseFromSourceMappedId')
        )
        .drop('efo_id')

        # Adding constants:
        .withColumn('datasourceId', lit('crispr'))
        .withColumn('datatypeId', lit('affected_pathway'))
        .persist()
    )

    logging.info(f'Saving {pooled_evidence_df.count()} CRISPR evidence in JSON format, to: {out_file}')

    (
        pooled_evidence_df
        .coalesce(1)
        .write.format('json').mode('overwrite').option('compression', 'gzip').save(out_file)
    )


if __name__ == "__main__":

    # Parse CLI arguments:
    parser = argparse.ArgumentParser(
        description='Parse essential cancer genes identified using CRISPR assays and prioritised in Project Score'
    )
    parser.add_argument('-d', '--descriptions_file',
                        help='Name of tsv file with the description of the method per cancer type',
                        type=str, required=True)
    parser.add_argument('-e', '--evidence_file',
                        help='Name of tsv file with the priority score',
                        type=str, required=True)
    parser.add_argument('-c', '--cell_types_file',
                        help='Name of tsv file with cell line names per cancer type',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file. (gzip compressed json)',
                        type=str, required=True)
    parser.add_argument('-l', '--log_file',
                        help='Name of log file. If not provided logs are written to standard error.',
                        type=str, required=False)

    args = parser.parse_args()
    desc_file = args.descriptions_file
    evid_file = args.evidence_file
    cell_file = args.cell_types_file
    out_file = args.output_file

    # If no logfile is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    main(desc_file, evid_file, cell_file, out_file)
