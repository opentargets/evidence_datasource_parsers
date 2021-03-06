#!/usr/bin/env python


'''
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=$HOME/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH

# Example input:
in_l2g = 'gs://genetics-portal-data/l2g/200127' # table with locus 2 gene scores
in_toploci = 'gs://genetics-portal-data/v2d/200207/toploci.parquet' # table with the top loci
in_study = 'gs://genetics-portal-data/v2d/200207/studies.parquet' # table with the studies
in_varindex = 'gs://genetics-portal-data/variant-annotation/190129/variant-annotation.parquet' # Table with the variants (gnomad)
in_csq_eco = 'gs://genetics-portal-data/lut/vep_consequences.tsv' # Table with consequence codes # (eco and label)

out_file = 'gs://genetics-portal-analysis/l2g-platform-export/data/l2g_joined.2020.02.28_exploded.parquet' # (output file.)

# Example call:
python GeneticsPortal_prepare_data.py --locus2gene gs://genetics-portal-data/l2g/200127 \
    --toploci gs://genetics-portal-data/v2d/200207/toploci.parquet \
    --study gs://genetics-portal-data/v2d/200207/studies.parquet \
    --variantIndex gs://genetics-portal-data/variant-annotation/190129/variant-annotation.parquet \
    --ecoCodes gs://genetics-portal-data/lut/vep_consequences.tsv' # Table with consequence codes \
    --output gs://genetics-portal-analysis/l2g-platform-export/data/l2g_joined.2020.02.28_exploded.parquet

'''

from datetime import datetime
import argparse
import os
import sys
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *

def load_eco_dict(inf):
    '''
    Loads the csq to eco scores into a dict
    Returns: dict
    '''

    # Load
    eco_df = (
        spark.read.csv(inf, sep='\t', header=True, inferSchema=True)
        .select('Term','Accession',col('eco_score').cast(DoubleType()))
    )

    # Convert to python dict
    eco_dict = {}
    eco_link_dict = {}
    for row in eco_df.collect():
        eco_dict[row.Term] = row.eco_score
        eco_link_dict[row.Term] = row.Accession
    
    return (eco_dict,eco_link_dict)

def main():

    ##
    ## Parsing parameters:
    ##
    parser = argparse.ArgumentParser(description='This script pulls together data from Open Targets Genetics portal to generate Platform evidences.')

    # Database related input:
    parser.add_argument('--locus2gene', help='Input table containing locus to gene scores.', type=str, required=True)
    parser.add_argument('--toploci', help='Table containing top loci for all studies.', type=str, required=True)
    parser.add_argument('--study', help='Table with all the studies.', type=str, required=True)
    parser.add_argument('--variantIndex', help='Table with the variant indices (from gnomad).', type=str, required=True)
    parser.add_argument('--ecoCodes', help='Table with consequence ECO codes.', type=str, required=True)
    parser.add_argument('--output', help='Output file name.', type=str, required=True)
    args = parser.parse_args()

    # Parse input parameters:  
    in_l2g = args.locus2gene
    in_toploci = args.toploci
    in_study = args.study
    in_varindex = args.variantIndex
    in_csq_eco = args.ecoCodes

    # Parse output parameter
    out_file = args.output
    out_path = '/'.join(out_file.split('/')[:-1]) # Will be used for logging.

    ##
    ## Initialize spark session
    ##
    global spark
    spark = (pyspark.sql.SparkSession.builder.getOrCreate())
    print('Spark version: ', spark.version)

    ##
    ## Load datasets
    ##

    # Load locus-to-gene (L2G) score data
    l2g = (
        spark.read.parquet(in_l2g)
        # Keep results trained on high or medium confidence gold-standards
        .filter(col('training_gs') == 'high_medium')
        # Keep results from xgboost model
        .filter(col('training_clf') == 'xgboost')
        # Only keep study, variant, gene and score info
        .select(
            'study_id',
            'chrom', 'pos', 'ref', 'alt',
            'gene_id',
            'y_proba_full_model',
        )
    )

    # Load association statistics (only pvalue is required) from top loci table
    pvals = (
        spark.read.parquet(in_toploci)
        # # Calculate pvalue from the mantissa and exponent
        # .withColumn('pval', col('pval_mantissa') * pow(10, col('pval_exponent')))
        # # NB. be careful as very small floats will be set to 0, we can se these
        # # to the smallest possible float instead
        # .withColumn('pval',
        #     when(col('pval') == 0, sys.float_info.min)
        #     .otherwise(col('pval'))
        # )
        # Keep required fields
        .select('study_id', 'chrom', 'pos', 'ref', 'alt', 
            'pval_mantissa', 'pval_exponent','odds_ratio','oddsr_ci_lower', 'oddsr_ci_upper')
    )

    # Load (a) disease information, (b) sample size from the study table
    study_info = (
        spark.read.parquet(in_study)
        .select(
            'study_id', 'pmid', 'pub_date', 'pub_author', 'trait_reported',
            'trait_efos',
            col('n_initial').alias('sample_size') # Rename to sample size
        )
                
        # Warning! Not all studies have an EFO annotated. Also, some have
        # multiple EFOs! We need to decide a strategy to deal with these.
        
        # # For example, only keep studies with 1 efo:
        # .filter(size(col('trait_efos')) == 1)
        # .withColumn('efo', col('trait_efos').getItem(0))
        # .drop('trait_efos')

        # Or, drop rows with no EFO and then explode array to multiple rows
        .withColumn('trait_efos', when(col('trait_efos').isNotNull(),
                                       expr('filter(trait_efos, t -> length(t) > 0)')))
        .filter(col('trait_efos').isNotNull() & (size(col('trait_efos')) > 0))
        .withColumn('efo', explode(col('trait_efos')))
        .drop('trait_efos')
    )

    # Get mapping for rsIDs:
    rsID_map = (
        spark.read.parquet(in_varindex)
        # chrom_b38|pos_b38
        # Explode consequences, only keeping canonical transcript
        .selectExpr(
            'chrom_b38 as chrom', 'pos_b38 as pos', 'ref', 'alt', 'rsid'
        )
    )

    # Load consequences:
    var_consequences = (
        spark.read.parquet(in_varindex)
        # chrom_b38|pos_b38
        # Explode consequences, only keeping canonical transcript
        .selectExpr(
            'chrom_b38 as chrom', 'pos_b38 as pos', 'ref', 'alt',
            'vep.most_severe_consequence as most_severe_csq',
            '''explode(
                filter(vep.transcript_consequences, x -> x.canonical == 1)
            ) as tc
            '''
        )
        # Keep required fields from consequences struct
        .selectExpr(
            'chrom', 'pos', 'ref', 'alt', 'most_severe_csq',
            'tc.gene_id as gene_id', 
            'tc.consequence_terms as csq_arr',
        )
    )

    ##
    ## Get most severe consequences:
    ##
    
    # Load term to eco score dict
    # (eco_dict,eco_link_dict) = spark.sparkContext.broadcast(load_eco_dict(in_csq_eco))
    eco_dicts = spark.sparkContext.broadcast(load_eco_dict(in_csq_eco))
    
    get_link = udf(
        lambda x: eco_dicts.value[1][x],
        StringType()
    )

    # Extract most sereve csq per gene.
    # Create UDF that reverse sorts csq terms using eco score dict, then select
    # the first item. Then apply UDF to all rows in the data.
    get_most_severe = udf(
        lambda arr: sorted(arr, key=lambda x: eco_dicts.value[0].get(x, 0), reverse=True)[0],
        StringType()
    )
    
    var_consequences = (
        var_consequences.withColumn('most_severe_gene_csq', get_most_severe(col('csq_arr')))
        .withColumn('consequence_link',get_link(col('most_severe_gene_csq')))
    )

    ##
    ## Join datasets together
    ##

    processed = (
        l2g
        # Join L2G to pvals, using study and variant info as key
        .join(pvals, on=['study_id', 'chrom', 'pos', 'ref', 'alt'])
        # Join this to the study info, using study_id as key
        .join(study_info, on='study_id')
        # Join transcript consequences:
        .join(var_consequences, on = ['chrom', 'pos', 'ref', 'alt', 'gene_id'], how = 'left')
        # Join rsIDs:
        .join(rsID_map, on = ['chrom','pos', 'ref', 'alt'], how = 'left')
        # Filling with missing values:
        .fillna({'most_severe_gene_csq' : 'intergenic_variant', 'consequence_link' : 'http://purl.obolibrary.org/obo/SO_0001628'})
    )

    # Write output
    (
        processed
        .write.parquet(
            out_file,
            mode='overwrite'
        )
    )

    return 0


if __name__ == '__main__':

    main()
