#!/usr/bin/env python
"""This module adds a more granular description of the phenotype observed in the ClinPGX evidence."""

import argparse
import json
import logging
import sys
from typing import List, Optional

from openai import OpenAI
import pyspark.sql.functions as f
from pyspark import SparkFiles
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StructType, StructField, StringType, ArrayType

from common.evidence import initialize_sparksession, write_evidence_strings
from common.ontology import add_efo_mapping

def parse_phenotype_with_gpt(genotype_text: str, openai_client: OpenAI, gpt_model: str = "gpt-3.5-turbo-1106") -> Optional[List[str]]:
    """Query the OpenAI API to extract the phenotype from the genotype text."""
    prompt = f"""
        Context: We want to analyze ClinPGx clinical annotations. Their data includes a column, "genotypeAnnotationText", which typically informs about efficacy,side effects, or patient response variability given a specific genotype. The data is presented in a lengthy and complex format, making it challenging to quickly grasp the key phenotypic outcomes.

        Aim: To parse the observed effect in a short string so that the effect can be easily interpreted at a glance. The goal is to extract the essence of the pharmacogenetic relationship. This extraction helps in summarizing the data for faster and more efficient analysis.

        Please analyze the following examples from the "genotypeAnnotationText" column and extract the key phenotype as a concise description. Format the result as a JSON array. Each JSON must only contain one field: "gptExtractedPhenotype".

        Examples for extraction:
        1. "Patients with the CTT/del genotype (one copy of the CFTR F508del variant) and cystic fibrosis may have increased response when treated with ivacaftor/tezacaftor combination as compared to patients with the CTT/CTT genotype." -> Expected extraction: "increased response"
        2. "Patients with the AC genotype may have increased risk for gastrointestinal toxicity with taxane and platinum regimens as compared to patients with the CC genotype." -> Expected extraction: "risk of gastrointestinal toxicity"
        3. "Patients with the rs2032582 AA genotype may be more likely to respond to tramadol treatment as compared to patients with the CC genotype." -> Expected extraction: "increased response"
        4. "Patients receiving methotrexate to treat acute lymphoblastic leukemia (ALL), and the rs4149056 TT genotype may be less likely to require glucarpidase treatment as compared to patients with the CC or CT genotypes." -> Expected extraction: "less likely to require glucarpidase"
        5. "Patients with the TT genotype and hormone insensitive breast cancer may experience increased risk of chemotherapy-induced amenorrhea when treated with goserelin or combinations of cyclophosphamide, docetaxel, doxorubicin, epirubicin, and fluorouracil compared to patients with the CT genotype." -> Expected extraction: "risk of chemotherapy-induced amenorrhea"
        6. "Patients with the GG genotype and cancer may have an increased risk for drug toxicity and an increased response to treatment with cisplatin or carboplatin as compared to patients with the AA or AG genotype. Other genetic and clinical factors may also influence a patient's risk for toxicity and response to platinum-based chemotherapy." -> Expected extraction: "drug toxicity" and "increased response"

        Based on these examples, please extract the phenotype from the following text:

        "{genotype_text}"
    """
    completion = openai_client.chat.completions.create(
        model=gpt_model,
        response_format={ "type": "json_object" },
        messages=[
            {"role": "system", "content": "you are an expert in clinical pharmacology designed to output JSON."},
            {"role": "user", "content": prompt},
        ],
        seed=42,
    )
    try:
        generated_text = completion.choices[0].message.content
        json_obj = json.loads(generated_text)
        return json_obj.get('gptExtractedPhenotype', [])
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def parse_phenotypes(texts_to_parse: List[str], openai_client: OpenAI) -> DataFrame:
    """Parse the phenotypes from the given texts by calling the OpenAI API."""
    results_dict = {}
    for text in texts_to_parse:
        result = parse_phenotype_with_gpt(text, openai_client)
        if isinstance(result, list):
            results_dict[text] = result
        elif isinstance(result, str):
            results_dict[text] = [result]
    return spark.createDataFrame(
        list(results_dict.items()),
        StructType(
            [
                StructField("genotypeAnnotationText", StringType(), True),
                StructField("phenotypeText", ArrayType(StringType()), True),
            ]
        ),
    )

def update_phenotypes_lut(
    new_phenotypes_df: DataFrame,
    extracted_phenotypes_df: DataFrame,
) -> DataFrame:
    """Adds the new phenotypes to the extracted phenotypes table."""
    return (
        extracted_phenotypes_df.unionByName(new_phenotypes_df)
    )


def generate_evidence(
    spark: SparkSession,
    pgx_evidence_df: DataFrame,
    extracted_phenotypes_df: DataFrame,
    cache_dir: str
) -> DataFrame:
    """This module overwrites the `phenotypeText` and adds `diseaseFromSourceMappedId` field in the PGx evidence dataset.
    
    The original `phenotypeText` comes from PGx directly, however it is tipically of little value for the user (more context in https://github.com/opentargets/curation/blob/master/docs/pharmacogenetics.md).

    Args:
        pgx_evidence_df: Dataframe with the PGx evidence submitted by EVA.
        extracted_phenotypes_df: Dataframe containing the phenotypes extracted from `genotypeAnnotationText`
        cache_dir: Directory to store the OnToma cache files in.
    """
    enriched_pgx_df = (
        pgx_evidence_df.drop("phenotypeText", "phenotypeFromSourceId").join(extracted_phenotypes_df, on="genotypeAnnotationText", how="left")
        .withColumn("phenotypeText", f.explode_outer("phenotypeText"))
        .distinct()
        .persist()
    )
    enriched_pgx_df = add_efo_mapping(
        enriched_pgx_df.select("*", f.col("phenotypeText").alias("diseaseFromSource"), f.lit(None).alias("diseaseFromSourceId")),
        spark,
        cache_dir
    ).withColumnRenamed("diseaseFromSourceMappedId", "phenotypeFromSourceId").drop("diseaseFromSource")
    logging.info("Disease mappings have been added.")

    assert enriched_pgx_df.select("studyId").distinct().count() >= pgx_df.select("studyId").distinct().count(), "Fewer ClinPGx references after processing."
    return enriched_pgx_df

def add_variantid_column(input_df: DataFrame) -> DataFrame:
    """Based on the content of the genotypeId column, adds a variantId column to the dataset"""
    return (
        input_df
        # split genotypeId column into chr pos ref alt columns
        .select("genotypeId", f.from_csv(f.col("genotypeId"), "chr string, pos string, ref string, alt string", {'sep': '_'}).alias("genotype_split"))
        .select("genotypeId", "genotype_split.*").toDF("genotypeId", "chr", "pos", "ref", "alt")
        # split alt column and explode
        .withColumn("alt", f.explode(f.split(f.col("alt"), ',')))
        .filter(~(f.col("ref") == f.col("alt")))
        .select("genotypeId", f.concat_ws('_', f.col("chr"), f.col("pos"), f.col("ref"), f.col("alt")).alias("variantId"))
        .join(input_df, on="genotypeId", how="right")
    )

def get_parser():
    """Get parser object for script ClinPGx.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--evidence_path",
        help="Input gzipped JSON with the PharmGKB evidence submitted by EVA",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--extracted_phenotypes_path",
        help="Input remote JSON file containing the phenotypes that have already been regenerated.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--openai_api_key", help="OpenAI API key. From https://platform.openai.com/api-keys.",
        type=str,
        required=False
    )
    parser.add_argument(
        "--output_evidence_path", help="Output gzipped json file containing the pharmgkb evidence with a new `phenotypeText` and `diseaseFromSourceMappedId` fields.",
        type=str,
        required=True
    )
    parser.add_argument(
        "--output_phenotypes_path", help="Output json file containing an updated version of the pharmgkb genotype to phenotype automatic annotation. The `curation` repo should be updated with this file.",
        type=str,
        required=False
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )
    return parser


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.StreamHandler(sys.stderr)

    spark = initialize_sparksession()
    args = get_parser().parse_args()

    # Read data
    logging.info(f"PGx evidence JSON file: {args.evidence_path}")
    logging.info(f"Table of genotype descriptions to extract traits: {args.extracted_phenotypes_path}")
    spark.sparkContext.addFile(args.extracted_phenotypes_path)
    pgx_phenotypes_df = spark.read.json(SparkFiles.get(args.extracted_phenotypes_path.split("/")[-1]))
    pgx_df = spark.read.json(args.pharmgkb_evidence_path)

    pgx_df = generate_evidence(
        spark,
        pgx_df,
        pgx_phenotypes_df,
        args.cache_dir
    )

    unparsed_texts = pgx_df.filter(f.col("phenotypeText").isNull()).select("genotypeAnnotationText").distinct()
    if unparsed_texts.count() == 0:
        logging.info("All evidence have been correctly parsed.")
    else:
        logging.error(f"There are  {unparsed_texts.count()} evidence without a phenotype. Please update the extracted phenotypes table before evidence generation.")
        client = OpenAI(api_key=args.openai_api_key)

        new_phenotypes_df = parse_phenotypes(
            unparsed_texts.toPandas()["genotypeAnnotationText"].to_list(),
            client
        )
        updated_phenotypes_df = update_phenotypes_lut(spark, new_phenotypes_df, pgx_phenotypes_df)
        updated_phenotypes_df.toPandas().to_json(args.output_phenotypes_path, orient="records")  # type: ignore
        logging.info(f"Updated phenotypes have been saved to {args.output_phenotypes_path}. Exiting.")
        pgx_df = generate_evidence(
            spark,
            pgx_df,
            updated_phenotypes_df,
            args.cache_dir
        )

    pgx_df = add_variantid_column(pgx_df)
    logging.info("Added variantId column.")
    write_evidence_strings(pgx_df, args.output_evidence_path)
    logging.info(f"{pgx_df.count()} evidence strings have been saved to {args.output_evidence_path}. Exiting.")

        

