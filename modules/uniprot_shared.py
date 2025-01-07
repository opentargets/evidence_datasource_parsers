"""Shared functions for the UniProt module."""

from __future__ import annotations

import logging
from typing import Any

from pyspark.sql import Column
from pyspark.sql import functions as f
from SPARQLWrapper import JSON, SPARQLWrapper

logger = logging.getLogger(__name__)


class UniprotShared:
    """Shared functions for the UniProt module."""

    UNIPROT_SPARQL_ENDPOINT = "https://sparql.uniprot.org/sparql"
    UNIPROT_SPARQL_QUERY = ""
    LOW_CONFIDENCE_FLAGS = [
        "The protein represented in this entry may be involved in disease pathogenesis",
        "Disease susceptibility may be associated with variants affecting the gene represented in this entry",
        "The disease may be caused by variants affecting distinct genetic loci, including the gene represented in this entry",
        "The disease may be caused by variants affecting the gene represented in this entry",
        "The gene represented in this entry may be involved in disease pathogenesis",
        "The gene represented in this entry may act as a disease modifier",
    ]
    SPARK_SESSION = None

    def __init__(self: UniprotShared, spark: SparkSession) -> None:
        """Initializes the class."""
        self.SPARK_SESSION = spark

    @classmethod
    def extract_uniprot_data(
        cls: type[UniprotShared],
    ) -> list[dict[str, Any]]:
        """Extracts Uniprot variants data from the Uniprot SPARQL API.

        Returns:
            list[dict[str, Any]]: The extracted Uniprot variants data.
        """
        logger.info(
            f"Extracting Uniprot variants data from {cls.UNIPROT_SPARQL_ENDPOINT} API endpoint."
        )
        # Initialize the SPARQL endpoint
        sparql = SPARQLWrapper(cls.UNIPROT_SPARQL_ENDPOINT)

        # Set the query and return format
        sparql.setQuery(cls.UNIPROT_SPARQL_QUERY)
        sparql.setReturnFormat(JSON)

        # Execute the query and fetch results
        results = sparql.query().convert()
        logger.info("Data extraction completed.")

        # Handling unexpected results:
        if not isinstance(results, dict):
            logger.error("No results found.")
            raise ValueError("SPARQL data extraction failed.")

        evidence_data = []
        for result in results["results"]["bindings"]:
            evidence = {}
            for key, value in result.items():
                # Extracting only MIM disease cross-references:
                if key == "disease_crossrefs" and "mim" not in value["value"]:
                    continue

                if value["type"] == "literal":
                    evidence[key] = value["value"]
                elif "mim" in value["value"]:
                    evidence[key] = value["value"]
                elif value["type"] == "uri":
                    evidence[key] = cls.get_uri_leaf(value["value"])

            evidence_data.append(evidence)

        logger.info(f"Number of disease/target/variant evidence: {len(evidence_data)}.")
        return evidence_data

    @staticmethod
    def get_uri_leaf(uri: str) -> str:
        """Extract the leaf of a URI.

        Args:
            uri (str): The URI.

        Returns:
            str: The leaf of the URI.
        """
        return uri.split("/")[-1]

    @classmethod
    def map_confidence(cls: type[UniprotShared], confidence_column: Column) -> Column:
        """Map confidence flags to a boolean column.

        Args:
            confidence_column (Column): The confidence column.

        Returns:
            Column: The mapped confidence column.
        """
        return (
            f.when(
                f.regexp_replace(f.split(confidence_column, r"\. ")[0], r"\.", "").isin(
                    cls.LOW_CONFIDENCE_FLAGS
                ),
                f.lit("low-confidence"),
            )
            .otherwise(f.lit("high-confidence"))
            .alias("confidence")
        )
