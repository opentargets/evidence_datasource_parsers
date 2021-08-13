#!/usr/bin/env python3
"""Evidence parser for the Cancer Biomarkers database."""

import argparse
import logging
import os
import tempfile

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import (
    array_distinct, coalesce, col, collect_set, concat, explode, expr, flatten, lit,
    initcap, regexp_extract, regexp_replace, split, struct, translate, trim, upper, when
)

ALTERATIONTYPE2FUNCTIONCSQ = {
    # TODO: Map BIA
    'MUT': 'SO_0001777', # somatic_variant
    'CNA': 'SO_0001563', # copy_number_change
    'FUS': 'SO_0001882', # feature_fusion
    'EXPR': 'SO_0001540', # level_of_transcript_variant
    'BIA': None
}

DRUGRESPONSE2EFO = {
    # TODO: Map Increased Toxicity and Resistant 
    'Responsive': 'GO_0042493', # response to drug
    'Not Responsive': 'HP_0020174', # Refractory drug response
    'Resistant': None,
    'Increased Toxicity': None,
    'Increased Toxicity (Myelosupression)': 'EFO_0007053', # myelosuppression
    'Increased Toxicity (Ototoxicity)': 'EFO_0006951', # ototoxicity
    'Increased Toxicity (Hyperbilirubinemia)': 'HP_0002904', # Hyperbilirubinemia
    'Increased Toxicity (Haemolytic Anemia)': 'EFO_0005558' # hemolytic anemia
}

class cancerBiomarkersEvidenceGenerator():
    
    def __init__(self):
        # Create spark session
        self.spark = (
            SparkSession.builder
            .appName('cancer-biomarkers')
            .getOrCreate()
        )
    
    def main(
        self,
        biomarkers_table: str,
        source_table: str,
        disease_table: str,
        drug_index: str,
        output_file: str
    ) -> None:
        """Ingests the processed tables to compute the evidence generation"""

    # Import data
        biomarkers_df = self.spark.read.csv(biomarkers_table, sep='\t', header=True)
        source_df = self.spark.read.json(source_table).select(
            col('label').alias('niceName'),
            'source', 'url')
        disease_df = self.spark.read.json(disease_table).select(
            regexp_replace(col('name'), '_', '').alias('tumor_type'),
            regexp_extract(col('url'), '[^/]+$', 0).alias('diseaseFromSourceMappedId'))
        drugs_df = self.spark.read.parquet(drug_index).select(
            col('id').alias('drugId'), col('name').alias('drug'))
        
        # Process tables
        biomarkers_enriched = self.process_biomarkers(biomarkers_df, source_df, disease_df, drugs_df)
        evidence = (biomarkers_enriched
            .withColumn('datasourceId', lit('cancer_genome_interpreter'))
            .withColumn('datatypeId', lit('affected_pathway'))
            # biomarker populated above
            .withColumnRenamed('tumor_type_full_name', 'diseaseFromSource')
            .withColumnRenamed('drug', 'drugFromSource')
            # diseaseFromSourceMappedId, drugId populated above
            .withColumnRenamed('Association', 'drugResponse')
            .withColumnRenamed('EvidenceLevel', 'confidence')
            .withColumnRenamed('alteration_type', 'variantFunctionalConsequenceId')
            # literature and urls populated above
            .withColumnRenamed('gene', 'targetFromSourceId')
            .withColumnRenamed('info', 'variantAminoacidDescriptions')
            .drop('tumor_type', 'source', 'Alteration', 'DrugFullName', 'niceName', 'url')
        )

        # Group evidence
        out_evidence = (evidence
            .groupBy('datasourceId', 'datatypeId', 'biomarker', 'drugFromSource', 
                    'drugResponse', 'targetFromSourceId', 'diseaseFromSource', 'diseaseFromSourceMappedId', 'drugId', 'confidence')
            .agg(
                collect_set('variantFunctionalConsequenceId').alias('variantFunctionalConsequenceId'),
                collect_set('literature').alias('literature'),
                collect_set('urls').alias('urls'),
                collect_set('variantAminoacidDescriptions').alias('variantAminoacidDescriptions')
            )
            .withColumn('variantAminoacidDescriptions', flatten(col('variantAminoacidDescriptions')))
            # Replace empty lists with null values
            .withColumn('literature', expr('filter(literature, x -> x is not null)'))
            .withColumn('urls', expr('filter(urls, x -> x is not null)'))
            .withColumn('variantAminoacidDescriptions', expr('filter(urls, x -> x is not null)'))
        )

        with tempfile.TemporaryDirectory() as tmp_dir_name:
                (
                    out_evidence.coalesce(1).write.format('json').mode('overwrite')
                    .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
                )
                json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
                assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
                os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)

    def process_biomarkers(
        self,
        biomarkers_df: DataFrame,
        source_df: DataFrame,
        disease_df: DataFrame,
        drugs_df: DataFrame
    ) -> DataFrame:
        """The diverse steps to prepare and enrich the input table are computed"""

        return (biomarkers_df
            .select(
                coalesce(col('IndividualMutation'), col('Biomarker')).alias('biomarker'),
                array_distinct(split(col('Gene'), ';')).alias('gene'),
                'Alteration',
                array_distinct(split(col('AlterationType'), ';')).alias('alteration_type'),
                array_distinct(split(col("PrimaryTumorTypeFullName"), ";")).alias('tumor_type_full_name'),
                array_distinct(split(col('Drug'), ';')).alias('drug'),
                'DrugFullName', 'Association', 'EvidenceLevel',
                array_distinct(split(col('Info'), ';')).alias('info'),
                array_distinct(split(col('Source'), ';')).alias('source')
            )
            .withColumn('tumor_type_full_name', explode(col('tumor_type_full_name')))
            .withColumn('tumor_type', translate(col("tumor_type_full_name"), " -", ""))
            .withColumn('gene', explode(col('gene')))
            .withColumn('drug', explode(col('drug')))
            # Override specific genes
            .withColumn('gene',
                  when(col('gene') == 'C15orf55', 'NUTM1')
                    .when(col('gene') == 'MLL', 'KMT2A')
                    .when(col('gene') == 'MLL2', 'KMT2D')
                    .otherwise(col('gene')))
            .withColumn('gene', upper(col('gene')))
            # Split source into literature and urls
            # literature contains PMIDs
            # urls are enriched from the source table if not a CT
            .withColumn('source', explode(col('source')))
            .withColumn('source', trim(regexp_extract(col('source'), '(PMID:\d+)|([\w ]+)', 0).alias('source')))
            .join(source_df, on='source', how='left')
            .withColumn('literature', when(col('source').startswith('PMID'), regexp_extract(col('source'), '(PMID:)(\d+)', 2)))
            .withColumn(
                'urls', 
                when(
                    col('source').startswith('NCT'),
                    struct(
                       lit('Clinical Trials').alias('niceName'),
                       concat(lit('https://clinicaltrials.gov/ct2/show/'), col('source')).alias('url')))
                .when(
                    (~col('source').startswith('PMID')) | (~col('source').startswith('NCIT')),
                    struct(col('niceName'), col('url'))))
            # Enrich data
            .withColumn('alteration_type', explode(col('alteration_type')))
            .replace(to_replace=ALTERATIONTYPE2FUNCTIONCSQ, subset=['alteration_type'])
            .replace(to_replace=DRUGRESPONSE2EFO, subset=['Association'])
            .join(disease_df, on='tumor_type', how='left')
            .withColumn('drug', upper(col('drug')))
            .withColumn(
                # drug class is coalesced when the precise name of the medicine is not provided 
                'drug',
                when(col('drug') == '[]', col('DrugFullName')).otherwise(col('drug')))
            .join(drugs_df, on='drug', how='left')
            .withColumn('drug', initcap(col('drug')))
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--biomarkers_table', help='Input TSV table containing association data for the biomarkers.', required=True)
    parser.add_argument('--source_table', help='Input JSON file with the annotations to reference the source of the associations.', required=True)
    parser.add_argument('--disease_table', help='Input JSON file with the mapping of the tumor types to EFO IDs.', required=True)
    parser.add_argument('--drug_index', help='Directory of parquet files that stores OT\'s disease index.', required=True)
    parser.add_argument('--output_file', help='Output gzipped json file containing the evidence.', required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    cancerBiomarkersEvidenceGenerator().main(
        args.biomarkers_table,
        args.source_table,
        args.disease_table,
        args.drug_index,
        args.output_file
    )