#!/usr/bin/env python
"""This script generates disease/target evidence based on INTOGEN data."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import TYPE_CHECKING

from pyspark.sql import DataFrame
from pyspark.sql import functions as f
from pyspark.sql import types as t

from common.evidence import (
    initialize_sparksession,
    write_evidence_strings,
)

if TYPE_CHECKING:
    from pyspark.sql import Column

DOI_TO_PMID_MAPPING = {
    "10.1038/ng.2529": "23334666",
    "10.1038/s41588-023-01321-1": "36928603",
    "10.1038/ng.3940": "28825729",
    "10.1101/162784": "32025007",
    "10.1172/JCI128227": "31483290",
    "10.1038/s41591-019-0582-4": "31570822",
    "10.1002/hep.27198": "24798001",
    "10.1038/ng.3520": "26928227",
    "10.1038/s41588-018-0078-z": "29610475",
}


def parse_source(source_column: Column) -> Column:
    processed = f.when(source_column.startswith("WEB"), f.lit(None)).otherwise(
        f.trim(f.split(source_column, ":")[1])
    )

    map_udf = f.udf(lambda v: DOI_TO_PMID_MAPPING.get(v, v), t.StringType())
    mapped = map_udf(processed)

    return f.when(mapped.isNotNull(), f.array(mapped))


class intogenEvidenceGenerator:
    def __init__(
        self: intogenEvidenceGenerator,
        inputGenes: str,
        inputCohorts: str,
        diseaseMapping: str | None,
    ) -> None:
        # Initialize spark session:
        self.spark = initialize_sparksession()

        # Reading cohort table:
        cohorts = self.spark.read.csv(
            inputCohorts, sep=r"\t", header=True, inferSchema=True
        ).drop("SAMPLES")

        # Reading cancer driver genes:
        self.evidence = (
            self.spark.read.csv(inputGenes, sep=r"\t", header=True, inferSchema=True)
            # Join with cohort table:
            .join(cohorts, on="COHORT", how="inner")
            .select(
                # Adding constant columns:
                f.lit("somatic_mutation").alias("datatypeId"),
                f.lit("intogen").alias("datasourceId"),
                # Extracting gene fields:
                f.col("SYMBOL").alias("targetFromSourceId"),
                # Extracting cohort specific annotation:
                f.col("COHORT").alias("cohortId"),
                f.col("COHORT_NICK").alias("cohortShortName"),
                f.col("COHORT_NAME").alias("cohortDescription"),
                # Splitting methods:
                f.split(f.col("METHODS"), ",").alias("significantDriverMethods"),
                # Generating mutated samples column:
                f.array(
                    f.struct(
                        # Parse functional consequences:
                        f.when(f.col("ROLE") == "Act", "SO_0002053")  # Gain of function
                        .when(f.col("ROLE") == "LoF", "SO_0002054")  # Loss of function
                        .otherwise(None)
                        .alias("functionalConsequenceId"),
                        # Extract samples:
                        f.col("SAMPLES").alias("numberSamplesTested"),
                        f.col("SAMPLES").alias("numberMutatedSamples"),
                    )
                ).alias("mutatedSamples"),
                # Adding disease name:
                f.col("CANCER").alias("Cancer_type_acronym"),
                f.col("CANCER_NAME").alias("diseaseFromSource"),
                # Adding score:
                f.col("QVALUE_COMBINATION").alias("resourceScore"),
                parse_source(f.col("REFERENCE")).alias("literature"),
            )
        )

        # Mapping step executed if mapping file is provided:
        if diseaseMapping is not None:
            logging.info(f"Applying disease mapping from {diseaseMapping}.")
            self.evidence = self.cancer2EFO(diseaseMapping)

            # Extracting stats about mapping:
            unmapped_evidence_count = self.evidence.filter(
                f.col("diseaseFromSourceMappedId").isNull()
            ).count()
            unmapped_disease_count = (
                self.evidence.filter(f.col("diseaseFromSourceMappedId").isNull())
                .select("diseaseFromSource")
                .distinct()
                .count()
            )
            logging.info(
                f"Number of evidence without EFO mapping: {unmapped_evidence_count}"
            )
            logging.info(
                f"Number of diseases without EFO mapping: {unmapped_disease_count}"
            )

        # Dropping cancer type acronym:
        self.evidence = self.evidence.drop("Cancer_type_acronym")

    def write_evidence(self: intogenEvidenceGenerator, evidence_file: str) -> None:
        write_evidence_strings(self.evidence, evidence_file)
        logging.info(
            f"{self.evidence.count()} evidence strings saved into {evidence_file}. Exiting."
        )

    def cancer2EFO(self: intogenEvidenceGenerator, diseaseMapping: str) -> DataFrame:
        diseaseMappingsFile = self.spark.read.csv(
            diseaseMapping, sep=r"\t", header=True
        ).select(
            "Cancer_type_acronym",
            f.trim(f.col("EFO_id")).alias("diseaseFromSourceMappedId"),
        )

        return self.evidence.join(
            diseaseMappingsFile, on="Cancer_type_acronym", how="left"
        )


def main(inputGenes, inputCohorts, diseaseMapping, outputFile):
    # Logging parameters
    logging.info(f"intOGen driver genes table: {inputGenes}")
    logging.info(f"intOGen cohorts table: {inputCohorts}")
    logging.info(f"Output file: {outputFile}")
    if diseaseMapping:
        logging.info(f"Cancer type to EFO ID table: {diseaseMapping}")
    else:
        logging.info("Disease mapping is skipped.")

    # Initialize evidence builder object
    evidenceBuilder = intogenEvidenceGenerator(inputGenes, inputCohorts, diseaseMapping)

    # Write evidence file:
    evidenceBuilder.write_evidence(outputFile)


if __name__ == "__main__":
    # Initiating parser
    parser = argparse.ArgumentParser(
        description="This script generates evidences for the intOGen data source."
    )

    parser.add_argument(
        "-g",
        "--inputGenes",
        required=True,
        type=str,
        help="Input source .tsv file listing the driver genes across the analyzed cohorts.",
    )
    parser.add_argument(
        "-c",
        "--inputCohorts",
        required=True,
        type=str,
        help="Input source .tsv file with information about the analyzed cohorts.",
    )
    parser.add_argument(
        "-d",
        "--diseaseMapping",
        required=False,
        type=str,
        help="Input look-up table containing the cancer type mappings to an EFO ID.",
    )
    parser.add_argument(
        "-o",
        "--outputFile",
        required=True,
        type=str,
        help="Gzipped JSON file containing the evidence strings.",
    )
    parser.add_argument(
        "-l",
        "--logFile",
        help="Destination of the logs generated by this script.",
        type=str,
        required=False,
    )

    # Parsing parameters
    args = parser.parse_args()
    inputGenes = args.inputGenes
    inputCohorts = args.inputCohorts
    diseaseMapping = args.diseaseMapping
    outputFile = args.outputFile

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format="%(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    main(inputGenes, inputCohorts, diseaseMapping, outputFile)
