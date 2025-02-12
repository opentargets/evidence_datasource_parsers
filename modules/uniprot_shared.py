"""Shared functions for the UniProt modules."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Any

from pyspark.sql import Column, DataFrame
from pyspark.sql import functions as f
from SPARQLWrapper import JSON, SPARQLWrapper

from common.ontology import add_efo_mapping

logger = logging.getLogger(__name__)


class UniprotShared(ABC):
    """Shared functions for the UniProt modules.

    Shared features between Uniprot variants and literature parsers:
    - Extract data from UniProt SPARQL API.
    - Map confidence flags to boolean columns.
    - Add EFO mappings to the UniProt evidence.
    - Extract unique leaf from a URI.
    - Get UniProt evidence - returns the evidence DataFrame.
    """

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
    EVIDENCE_COLUMNS: list[str] = []

    evidence_dataframe: DataFrame | None = None

    @abstractmethod
    def extract_evidence_from_uniprot(
        self: UniprotShared,
    ) -> UniprotShared:
        """Extracts evidence from UniProt.

        Returns:
            UniprotShared: The UniProt shared instance.
        """
        pass

    def extract_uniprot_data(
        self: UniprotShared,
    ) -> list[dict[str, Any]]:
        """Extracts Uniprot variants data from the Uniprot SPARQL API.

        Returns:
            list[dict[str, Any]]: The extracted Uniprot variants data.
        """
        logger.info(
            f"Extracting Uniprot variants data from {self.UNIPROT_SPARQL_ENDPOINT} API endpoint."
        )
        # Initialize the SPARQL endpoint
        sparql = SPARQLWrapper(self.UNIPROT_SPARQL_ENDPOINT)

        # Set the query and return format
        sparql.setQuery(self.UNIPROT_SPARQL_QUERY)
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
                    evidence[key] = self.get_uri_leaf(value["value"])

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
                f.lit("medium"),
            )
            .otherwise(f.lit("high"))
            .alias("confidence")
        )

    def get_evidence(self: UniprotShared, debug: bool = False) -> DataFrame:
        """Get the Uniprot evidence.

        Returning the expected columns for the evidence, however, if debug is set to True, all columns are returned.
        Serving as a debugging tool to explore the evidence data in detail.

        Args:
            debug (bool, optional): Whether to return all columns for debugging purposes. Defaults to False.

        Returns:
            DataFrame: The Uniprot evidence.

        Raises:
            ValueError: If no evidence data is found.
        """
        if not self.evidence_dataframe:
            raise ValueError("No evidence data found.")

        # Return all columns for debugging purposes:
        select_expression = "*" if debug else self.EVIDENCE_COLUMNS

        return self.evidence_dataframe.select(*select_expression)

    def add_efo_mapping(self: UniprotShared, ontoma_cache_dir: str) -> UniprotShared:
        """Add EFO mappings to the Uniprot evidence.

        Args:
            ontoma_cache_dir (str): The OnToma cache directory.

        Returns:
            UniprotVariantsParser: The UniprotVariantsParser.

        Raises:
            ValueError: If no evidence data is found or Spark session not initialized.
        """
        if not self.evidence_dataframe or not self.SPARK_SESSION:
            raise ValueError("No evidence data found or Spark session not initialized.")

        # Add EFO mappings:
        logger.info("Adding EFO mappings.")
        self.evidence_dataframe = add_efo_mapping(
            self.evidence_dataframe, self.SPARK_SESSION, ontoma_cache_dir
        )

        return self

    @staticmethod
    def extract_omim_crossref(
        omim_refs: Column,
    ) -> Column:
        """Extract OMIM cross-references from the OMIM references column.

        Args:
            omim_refs (Column): The OMIM references column.

        Returns:
            Column: The OMIM cross-references.
        """
        return f.concat(
            f.lit("OMIM:"),
            f.split(omim_refs, "/")[f.size(f.split(omim_refs, "/")) - 1],
        )

    @staticmethod
    def extract_pmids(literature: Column) -> Column:
        """Extract PMIDs from the literature column.

        Args:
            literature (Column): The literature column.

        Returns:
            Column: The PMIDs.
        """
        return f.transform(
            # Literature references are separated by a comma:
            f.split(literature, ", "),
            # Each element needs to be further split, and extract the last element:
            lambda uri: f.split(uri, "/")[f.size(f.split(uri, "/")) - 1],
        )
