"""Module to generate SLAPenrich disease/target evidence."""

from __future__ import annotations

import argparse
import logging
from typing import TYPE_CHECKING

import pyspark.sql.functions as f

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)

if TYPE_CHECKING:
    from pyspark.sql import DataFrame

logger = logging.getLogger(__name__)


class SLAPEnrichEvidenceParser:
    """Logic to generate SLAPenrich evidence."""

    evidence_data: DataFrame | None = None

    def __init__(
        self: SLAPEnrichEvidenceParser,
        input_file: str,
        disease_mapping_file: str,
        p_value_threshold: float = 1e-4,
    ) -> None:
        """Initialise parser object.

        Args:
            input_file (str): SLAPenrich input file with evidence
            disease_mapping_file (str): file from disease mapping is read
            p_value_threshold (float): p-value threshold applied on evidence. Default: 1e-4
        """
        logger.info("SLAPenrich parser was called with the following parameters:")
        logger.info(f"Input file: {input_file}")
        logger.info(f"Disease mapping file: {disease_mapping_file}")
        logger.info(f"P-value threshold: {p_value_threshold}")

        spark = initialize_sparksession()

        # Reading raw evidence:
        raw_evidence: DataFrame = spark.read.csv(
            input_file, sep=r"\t", header=True, inferSchema=True
        )

        logger.info(f"Number of raw evidnece: {raw_evidence.count()}")

        # Read disease mapping file:
        disease_mapping: DataFrame = spark.read.csv(
            disease_mapping_file, sep=r"\t", header=True
        ).select(f.col("Cancer_type_acronym").alias("ctype"), "EFO_id")

        # Read evidence file:
        self.evidence_data = (
            raw_evidence
            # Drop sub-significant rows:
            .filter(f.col("SLAPEnrichPval") <= p_value_threshold)
            # Join with disease mappings:
            .join(disease_mapping, on="ctype", how="left")
            # Extract relevant columns:
            .select(
                f.lit("slapenrich").alias("datasourceId"),
                f.lit("affected_pathway").alias("datatypeId"),
                f.array(f.lit("29713020")).alias("literature"),
                f.col("ctype").alias("diseaseFromSource"),
                f.col("gene").alias("targetFromSourceId"),
                f.col("SLAPEnrichPval").alias("resourceScore"),
                f.col("EFO_id").alias("diseaseFromSourceMappedId"),
                f.array(
                    f.struct(
                        f.split(f.col("pathway"), ": ").getItem(0).alias("id"),
                        f.split(f.col("pathway"), ": ").getItem(1).alias("name"),
                    )
                ).alias("pathways"),
            )
        )

        logger.info(f"Number of parsed evidence: {self.evidence_data.count()}")

    def get_evidence(self):
        """Retrun evidence data.

        Returns:
            DataFrame: parsed evidence data
        """
        return self.evidence_data

    def save_evidence(self, output_file_name: str) -> None:
        """Save evidence as gzipped json.

        Args:
            output_file_name (str): the name of the output evidence file.
        """
        write_evidence_strings(self.evidence_data, output_file=output_file_name)
        logger.info("Saving evidence strings finished.")


def main(input_file: str, disease_mapping_file: str, output_file: str) -> None:
    # Initialiase parser object:
    parser_object = SLAPEnrichEvidenceParser(
        input_file=input_file, disease_mapping_file=disease_mapping_file
    )

    # Save data:
    parser_object.save_evidence(output_file_name=output_file)


if __name__ == "__main__":
    # Initiating parser
    parser = argparse.ArgumentParser(
        description="This script generates evidences for the SLAPEnrich data source."
    )

    parser.add_argument(
        "-i", "--input_file", required=True, type=str, help="Input source .tsv file."
    )
    parser.add_argument(
        "-d",
        "--disease_mapping",
        required=False,
        type=str,
        help="Input look-up table containing the cancer type mappings to an EFO ID.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        type=str,
        help="Gzipped JSON file containing the evidence strings.",
    )

    # Parsing parameters
    args = parser.parse_args()
    input_file = args.input_file
    disease_mapping = args.disease_mapping
    output_file = args.output_file

    # Initialize logging:
    initialize_logger(__name__)

    main(input_file, disease_mapping, output_file)
