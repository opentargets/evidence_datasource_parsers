#!/usr/bin/env python
"""This module adds a more granular description of the phenotype observed in the PharmGKB evidence."""

import argparse
import logging
import sys
import json

from openai import OpenAI
import pyspark.sql.functions as f
from pyspark import SparkFiles
from pyspark.sql import SparkSession

from common.evidence import initialize_sparksession, write_evidence_strings
from common.ontology import add_efo_mapping

def parse_phenotype(genotype_text: str, openai_client: OpenAI, gpt_model: str = "gpt-3.5-turbo-1106") -> list | None:
    prompt = f"""
        Context: We want to analyze PharmGKB clinical annotations. Their data includes a column, "genotypeAnnotationText", which typically informs about efficacy,side effects, or patient response variability given a specific genotype. The data is presented in a lengthy and complex format, making it challenging to quickly grasp the key phenotypic outcomes.

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

def main(spark: SparkSession, pharmgkb_evidence_path: str, extracted_phenotypes_path: str, output_file_path: str, cache_dir: str) -> None:
    """This module overwrites the `phenotypeText` and adds `diseaseFromSourceMappedId` field in the PharmGKB evidence dataset.
    
    The original `phenotypeText` comes from PharmGKB directly, however it is tipically of little value for the user (more context in https://github.com/opentargets/curation/blob/master/docs/pharmacogenetics.md).

    Args:
        pharmgkb_evidence_path: Input gzipped JSON with the evidence submitted by ChEMBL.
        extracted_phenotypes_path: Input JSON with the phenotypes extracted from `genotypeAnnotationText`
        output_file_path: Output gzipped json file containing the pharmgkb evidence with a new `phenotypeText` and `diseaseFromSourceMappedId` fields.
        cache_dir: Directory to store the OnToma cache files in.
    """
    # Extract
    logging.info(f"PharmGKB evidence JSON file: {pharmgkb_evidence_path}")
    logging.info(f"Table of genotype descriptions to extracted phenotypes: {extracted_phenotypes_path}")
    spark.sparkContext.addFile(extracted_phenotypes_path)
    pharmgkb_df = spark.read.json(pharmgkb_evidence_path)
    pgx_phenotypes_df = spark.read.json(SparkFiles.get(extracted_phenotypes_path.split("/")[-1]))


    # Transform
    enriched_pharmgkb_df = (
        pharmgkb_df.drop("phenotypeText", "phenotypeFromSourceId").join(pgx_phenotypes_df, on="genotypeAnnotationText", how="left")
        .withColumn("phenotypeText", f.explode("phenotypeText"))
        .persist()
    )
    assert (enriched_pharmgkb_df.filter(f.col("phenotypeText").isNull()).count() == 0), "There are evidence without a phenotype. Please update the extracted phenotypes table before evidence generation."
    enriched_pharmgkb_df = add_efo_mapping(
        enriched_pharmgkb_df.select("*", f.col("phenotypeText").alias("diseaseFromSource"), f.lit(None).alias("diseaseFromSourceId")),
        spark,
        cache_dir
    ).withColumnRenamed("diseaseFromSourceMappedId", "phenotypeFromSourceId").drop("diseaseFromSource")
    logging.info("Disease mappings have been added.")

    # Load
    assert enriched_pharmgkb_df.count() >= pharmgkb_df.count(), "There are fewer evidence after processing."
    write_evidence_strings(enriched_pharmgkb_df, output_file_path)
    logging.info(f"{enriched_pharmgkb_df.count()} evidence strings have been saved to {output_file_path}. Exiting.")


def get_parser():
    """Get parser object for script pharmgkb.py."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--pharmgkb_evidence_path",
        help="Input gzipped JSON with the PharmGKB evidence submitted by EVA",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--extracted_phenotypes_path",
        help="Input TSV containing the categories of the clinical trial reason to stop. Instructions for applying the ML model here: https://github.com/ireneisdoomed/stopReasons.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--openai_api_key", help="OpenAI API key. From https://platform.openai.com/api-keys.",
        type=str,
        required=True
    )
    parser.add_argument(
        "--output_file_path", help="Output gzipped json file following the target safety liabilities data model.",
        type=str,
        required=True
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

    main(
        spark,
        args.pharmgkb_evidence_path,
        args.extracted_phenotypes_path,
        args.openai_api_key,
        args.output_file_path,
        args.cache_dir
    )
