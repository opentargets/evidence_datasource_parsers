#!/usr/bin/env python3
"""Evidence parser for the Cancer Biomarkers database."""

import argparse
import logging
import os
import re
import tempfile

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.functions import (
    array_distinct, coalesce, col, collect_set, concat, element_at, explode, flatten, lit,
    initcap, regexp_extract, regexp_replace, size, split, struct, translate, trim, udf, upper, when
)
from pyspark.sql.types import StringType, ArrayType

ALTERATIONTYPE2FUNCTIONCSQ = {
    # TODO: Map BIA
    'MUT': 'SO_0001777',  # somatic_variant
    'CNA': 'SO_0001563',  # copy_number_change
    'FUS': 'SO_0001882',  # feature_fusion
    'EXPR': None,
    'BIA': None
}

DRUGRESPONSE2EFO = {
    # TODO: Map Increased Toxicity and Resistance
    'Responsive': 'GO_0042493',  # response to drug
    'Not Responsive': 'HP_0020174',  # Refractory drug response
    'Resistant': None,
    'Increased Toxicity': None,
    'Increased Toxicity (Myelosupression)': 'EFO_0007053',  # myelosuppression
    'Increased Toxicity (Ototoxicity)': 'EFO_0006951',  # ototoxicity
    'Increased Toxicity (Hyperbilirubinemia)': 'HP_0002904',  # Hyperbilirubinemia
    'Increased Toxicity (Haemolytic Anemia)': 'EFO_0005558'  # hemolytic anemia
}


class cancerBiomarkersEvidenceGenerator():

    def __init__(self):
        # Create spark session
        self.spark = (
            SparkSession.builder
            .appName('cancer-biomarkers')
            .getOrCreate()
        )

        self.evidence = None

        self.get_variantId_udf = udf(
            cancerBiomarkersEvidenceGenerator.get_variantId,
            StringType()
        )
        self.zip_alterations_with_type_udf = udf(
            cancerBiomarkersEvidenceGenerator.zip_alterations_with_type,
            ArrayType(ArrayType(StringType()))
        )

    def main(
        self,
        biomarkers_table: str,
        source_table: str,
        disease_table: str,
        drug_index: str,
        output_file: str
    ) -> None:
        """Loads and processes inputs to generate the Cancer Biomarkers evidence strings"""

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

        # Process inputs to generate evidence strings
        evidence = self.compute_biomarkers(
            biomarkers_df, source_df, disease_df, drugs_df
        )

        # Write evidence strings
        self.write_evidence_strings(output_file)
        logging.info(f'{evidence.count()} evidence strings have been saved to {output_file}.')

    def compute_biomarkers(
        self,
        biomarkers_df: DataFrame,
        source_df: DataFrame,
        disease_df: DataFrame,
        drugs_df: DataFrame
    ) -> DataFrame:
        """The diverse steps to prepare and enrich the input table are computed"""

        biomarkers_enriched = (
            biomarkers_df
            .select(
                'Biomarker', 'IndividualMutation',
                array_distinct(split(col('Alteration'), ';')).alias('alterations'),
                array_distinct(split(col('Gene'), ';')).alias('gene'),
                split(col('AlterationType'), ';').alias('alteration_types'),
                array_distinct(split(col("PrimaryTumorTypeFullName"), ";")).alias('tumor_type_full_name'),
                array_distinct(split(col('Drug'), ';|,')).alias('drug'),
                'DrugFullName', 'Association', 'gDNA',
                array_distinct(split(col('EvidenceLevel'), ',')).alias('confidence'),
                array_distinct(split(col('Source'), ';')).alias('source')
            )
            .withColumn('confidence', explode(col('confidence')))
            .withColumn('tumor_type_full_name', explode(col('tumor_type_full_name')))
            .withColumn('tumor_type', translate(col('tumor_type_full_name'), ' -', ''))
            .withColumn('gene', explode(col('gene')))
            .withColumn('drug', explode(col('drug')))
            .withColumn('drug', translate(col('drug'), '[]', ''))
            # Override specific genes
            .withColumn(
                'gene',
                when(col('gene') == 'C15orf55', 'NUTM1')
                .when(col('gene') == 'MLL', 'KMT2A')
                .when(col('gene') == 'MLL2', 'KMT2D')
                .otherwise(col('gene'))
            )
            .withColumn('gene', upper(col('gene')))
            # Disambiguate alteration_type when biomarker consists of multiple alterations
            .withColumn(
                'tmp',
                self.zip_alterations_with_type_udf(col('alterations'), col('alteration_types')))
            .withColumn('tmp', explode(col('tmp')))
            .withColumn('alteration_type', element_at(col('tmp'), 2))
            .withColumn(
                'alteration',
                when(
                    ~col('IndividualMutation').isNull(),
                    col('IndividualMutation')
                )
                .otherwise(element_at(col('tmp'), 1))
            )
            .drop('tmp')
            # Clean special cases on the alteration string
            .withColumn(
                'alteration',
                when(
                    col('alteration').contains(':.'),
                    translate(col('alteration'), ':.', '')
                )
                .when(
                    col('alteration') == 'NRAS:.12.,.13.,.59.,.61.,.117.,.146.',
                    col('Biomarker')  # 'NRAS (12,13,59,61,117,146)'
                )
                .when(
                    # Fusion genes are described with '__'
                    # biomarker is a cleaner representation when there's one alteration
                    (col('alteration').contains('__')) & (~col('Biomarker').contains('+')),
                    col('Biomarker')
                )
                .otherwise(col('alteration'))
            )
            # Split source into literature and urls
            # literature contains PMIDs
            # urls are enriched from the source table if not a CT
            .withColumn('source', explode(col('source')))
            .withColumn('source', trim(regexp_extract(col('source'), r'(PMID:\d+)|([\w ]+)', 0).alias('source')))
            .join(source_df, on='source', how='left')
            .withColumn(
                'literature',
                when(col('source').startswith('PMID'), regexp_extract(col('source'), r'(PMID:)(\d+)', 2))
            )
            .withColumn(
                'urls',
                when(
                    col('source').startswith('NCT'),
                    struct(
                        lit('Clinical Trials').alias('niceName'),
                        concat(lit('https://clinicaltrials.gov/ct2/show/'), col('source')).alias('url')
                    )
                )
                .when(
                    (~col('source').startswith('PMID')) | (~col('source').startswith('NCIT')),
                    struct(col('niceName'), col('url'))
                )
            )
            # The previous conditional clause creates a struct regardless of
            # whether any condition is met. The empty struct is replaced with null
            .withColumn('urls', when(~col('urls.niceName').isNull(), col('urls')))
            # Enrich data
            .withColumn('variantFunctionalConsequenceId', col('alteration_type'))
            .replace(to_replace=ALTERATIONTYPE2FUNCTIONCSQ, subset=['variantFunctionalConsequenceId'])
            .replace(to_replace=DRUGRESPONSE2EFO, subset=['Association'])
            .join(disease_df, on='tumor_type', how='left')
            .withColumn('drug', upper(col('drug')))
            .withColumn(
                # drug class is coalesced when the precise name of the medicine is not provided
                'drug',
                when(col('drug') == '[]', col('DrugFullName')).otherwise(col('drug')))
            .join(drugs_df, on='drug', how='left')
            .withColumn('drug', initcap(col('drug')))
            # Translate variantId
            .withColumn(
                'variantId',
                when(~col('gDNA').isNull(), self.get_variantId_udf(col('gDNA')))
            )
            # Assign a GO ID when a gene expression data is reported
            .withColumn(
                'geneExpressionId',
                when(
                    (col('alteration_type') == 'EXPR') & (col('alteration').contains('over')),
                    'GO_0010628'
                )
                .when(
                    (col('alteration_type') == 'EXPR') & (col('alteration').contains('under')),
                    'GO_0010629'
                )
                .when(
                    (col('alteration_type') == 'EXPR') & (col('alteration').contains('norm')),
                    'GO_0010467'
                )
            )
            # Create variant struct
            .withColumn(
                'variant',
                when(
                    col('alteration_type') != 'EXPR',
                    struct(
                        col('alteration').alias('name'),
                        col('variantId').alias('id'),
                        col('variantFunctionalConsequenceId').alias('functionalConsequenceId')
                    )
                )
            )
            # Create geneExpression struct
            .withColumn(
                'geneExpression',
                when(
                    col('alteration_type') == 'EXPR',
                    struct(
                        col('alteration').alias('name'),
                        col('geneExpressionId').alias('id'))
                )
            )
        )

        pre_evidence = (
            biomarkers_enriched
            .withColumn('datasourceId', lit('cancer_genome_interpreter'))
            .withColumn('datatypeId', lit('affected_pathway'))
            .withColumnRenamed('tumor_type_full_name', 'diseaseFromSource')
            .withColumnRenamed('drug', 'drugFromSource')
            # diseaseFromSourceMappedId, drugId populated above
            .withColumnRenamed('Association', 'drugResponse')
            # confidence, literature and urls populated above
            .withColumnRenamed('gene', 'targetFromSourceId')
            .withColumnRenamed('Biomarker', 'biomarkerName')
            # variant, geneExpression populated above
            .drop(
                'tumor_type', 'source', 'alteration', 'alteration_type', 'IndividualMutation', 'geneExpressionId',
                'gDNA', 'variantFunctionalConsequenceId', 'variantId', 'DrugFullName', 'niceName', 'url')
        )

        # Group evidence
        self.evidence = (
            pre_evidence
            .groupBy('datasourceId', 'datatypeId', 'drugFromSource', 'drugId',
                     'drugResponse', 'targetFromSourceId', 'diseaseFromSource', 
                     'diseaseFromSourceMappedId', 'confidence', 'biomarkerName')
            .agg(
                collect_set('literature').alias('literature'),
                collect_set('urls').alias('urls'),
                collect_set('variant').alias('variant'),
                collect_set('geneExpression').alias('geneExpression'),
            )
            # Replace empty lists with null values
            .withColumn('literature', when(size(col('literature')) == 0, lit(None)).otherwise(col('literature')))
            .withColumn('urls', when(size(col('urls')) == 0, lit(None)).otherwise(col('urls')))
            .withColumn('variant', when(size(col('variant')) == 0, lit(None)).otherwise(col('variant')))
            .withColumn(
                'geneExpression',
                when(size(col('geneExpression')) == 0, lit(None))
                .otherwise(col('geneExpression')))
            # Collect variant info into biomarkers struct
            .withColumn(
                'biomarkers',
                struct(
                    'variant',
                    'geneExpression'
                ))
            .drop('variant', 'geneExpression')
            .distinct()
        )

        return self.evidence
    
    def write_evidence_strings(self, output_file: str) -> None:
        '''Exports the table to a compressed JSON file containing the evidence strings'''
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            (
                self.evidence.coalesce(1).write.format('json').mode('overwrite')
                .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
            )
            json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
            assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
            os.rename(os.path.join(tmp_dir_name, json_chunks[0]), output_file)

    @staticmethod
    def get_variantId(gDNA: str) -> str:
        '''
        Converts the genomic coordinates to the CHROM_POS_REF_ALT notation
        Ex.: 'chr14:g.105243048G_T' --> '14_105243048_G_T'
        '''
        translate_dct = {'chr': '', ':g.': '_', '>': '_', 'del': '_', 'ins': '_'}
        try:
            for k, v in translate_dct.items():
                gDNA = gDNA.replace(k, v)
            x, head, tail = re.split('^(.*?_\d+)', gDNA)
            if bool(re.search('\d+', tail)):
                tail = re.split('^(_\d+_)', tail)[-1]
            return head + '_' + tail
        except AttributeError:
            return

    @staticmethod
    def zip_alterations_with_type(alterations, alteration_type):
        '''
        Zips in a tuple the combination of the alteration w/ its correspondent type 
        so that when multiple alterations are reported, these can be disambiguated.
        By expanding the array of alteration types it accounts for the cases when
        several alterations are reported but only one type is given
        '''
        alteration_types = alteration_type * len(alterations)
        return list(zip(alterations, alteration_types))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--biomarkers_table',
                        help='Input TSV table containing association data for the biomarkers.', required=True)
    parser.add_argument('--source_table', required=True,
                        help='Input JSON file with the annotations to reference the source of the associations.')
    parser.add_argument('--disease_table',
                        help='Input JSON file with the mapping of the tumor types to EFO IDs.', required=True)
    parser.add_argument('--drug_index', help='Directory of parquet files that stores OT\'s disease index.',
                        required=True)
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
