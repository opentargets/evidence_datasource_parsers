"""Parser to build disease/target evidence based on Uniprot literature data."""

from __future__ import annotations

import argparse
import logging
from typing import Any

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as f
from SPARQLWrapper import JSON, SPARQLWrapper

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping


class UniprotLiteratureExtractor:
    """Class to extract Uniprot literature annotation from Unirot SPARQL API."""

    UNIPROT_SPARQL_ENDPOINT = "https://sparql.uniprot.org/sparql"

    UNIPROT_SPARQL_QUERY = """
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT
            ?protein 
            ?db
            ?disease_label
            ?disease_comment
            (GROUP_CONCAT(DISTINCT ?source; SEPARATOR=", ") AS ?sources)
        
        WHERE
        {
            # Defining the universe of proteins:
            # BIND(<http://purl.uniprot.org/uniprot/B1AK53> AS ?protein) .
            ?protein a up:Protein ;
                up:organism taxon:9606 .

            ?protein rdfs:seeAlso ?db .
            ?db up:database <http://purl.uniprot.org/database/MIM> .
            ?db rdfs:comment ?text .
            FILTER(CONTAINS(str(?text), "phenotype")) .
            ?disease_db rdfs:seeAlso ?db .
            ?disease_db rdf:type up:Disease .
            ?disease_db skos:prefLabel ?disease_label .
            ?disease_db rdfs:comment ?disease_comment .
            OPTIONAL {
                ?linkToEvidence rdf:object ?disease_db ;
                                up:attribution ?attribution .
                ?attribution up:source ?source .
                ?source a up:Journal_Citation .
            }       
        }

        GROUP BY 
            ?protein 
            ?db
            ?disease_label
            ?disease_comment
    """

    @classmethod
    def extract_literature_data(
        cls: type[UniprotLiteratureExtractor],
    ) -> list[dict[str, Any]]:
        """Extract Uniprot literature data from the Uniprot SPARQL API.

        Args:
            limit (int, optional): The number of rows to limit the data to. Defaults to 1000.

        Returns:
            list[dict[str, Any]]: The Uniprot literature data.
        """
        # Initialize the SPARQL endpoint
        sparql = SPARQLWrapper(cls.UNIPROT_SPARQL_ENDPOINT)

        # Set the query and return format
        sparql.setQuery(cls.UNIPROT_SPARQL_QUERY)
        sparql.setReturnFormat(JSON)

        # Execute the query and fetch results
        results = sparql.query().convert()
        logger.info("Data extraction completed.")

        literature_data = []
        for result in results["results"]["bindings"]:
            disease = {}
            for key, value in result.items():
                # Extracting only MIM disease cross-references:
                if key == "df" and "mim" not in value["value"]:
                    continue

                if value["type"] == "literal":
                    disease[key] = value["value"]
                elif "mim" in value["value"]:
                    disease[key] = value["value"]
                elif value["type"] == "uri":
                    disease[key] = cls.get_uri_leaf(value["value"])
            literature_data.append(disease)

        logger.info(
            f"Number of disease/target/variant evidence: {len(literature_data)}."
        )
        return literature_data

    @classmethod
    def get_dataframe(cls, spark: SparkSession) -> DataFrame:
        """Convert VEP data to a Spark DataFrame and minimally process it.

        Args:
            spark (SparkSession): The Spark session.

        Returns:
            DataFrame: The processed VEP data.
        """
        literature_data = cls.extract_literature_data()

        # Convert to a Spark DataFrame:
        return spark.createDataFrame(literature_data).select(
            # Extract disease information:
            f.col("disease_label").alias("diseaseFromSource"),
            f.concat(
                f.lit("OMIM:"),
                f.split("db", "/").getItem(f.size(f.split("db", "/")) - 1),
            ).alias("diseaseFromSourceId"),
            # Extract protein information:
            f.col("protein").alias("targetFromSourceId"),
            # Extract literature information:
            f.transform(
                f.filter(
                    f.split(f.col("sources"), ", "),
                    lambda source: f.length(source) > 0,
                ),
                lambda uri: f.split(uri, "/").getItem(f.size(f.split(uri, "/")) - 1),
            ).alias("literature"),
            f.col("disease_comment").alias("diseaseFromSourceDescription"),
        )

    @staticmethod
    def get_uri_leaf(uri: str) -> str:
        return uri.split("/")[-1]


def main(
    output_file: str,
    ontoma_cache_dir: str,
) -> None:
    # Get logger:
    logger.info("Starting Uniprot evidence parser.")
    logger.info(f"Output file: {output_file}")

    # Initialise spark session:
    spark = initialize_sparksession()

    # Extracting Uniprot evidence:
    uniprot_variants = UniprotLiteratureExtractor.get_dataframe(spark)

    # Add EFO mappings:
    logger.info("Adding EFO mappings.")
    evidence_df = add_efo_mapping(uniprot_variants, spark, ontoma_cache_dir)

    # Write data:
    logger.info("Writing data.")
    evidence_df.write.mode("overwrite").parquet("uniprot_debug.parquet")
    write_evidence_strings(
        (
            evidence_df
            # TODO: Based on targetToDiseaseAnnotation classify confidence
            # TODO: Based on variantToTargetAnnotation classify direction of effect/target modulation
            # The reference and altenate amino acids are there for troubleshooting purposes.
            .drop(
                "diseaseFromSourceDescription",
            ).withColumns(
                {
                    "datasourceId": f.lit("uniprot_literature"),
                    "datatypeId": f.lit("genetic_literature"),
                    "confidence": f.lit("high"),
                    "targetModulation": f.lit("up_or_down"),
                }
            )
        ),
        output_file,
    )


def parse_command_line_arguments() -> argparse.Namespace:
    """Parsing command line arguments for the Uniprot evidence parser.

    Returns:
        argparse.Namespace: The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Parser to build disease/target evidence based on Uniprot literature data."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="Path to the output file.",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        required=False,
        help="Path to the log file.",
    )
    parser.add_argument(
        "--ontoma_cache_dir",
        type=str,
        help="Path to the OnToma cache directory.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    # Exract command line arguments:
    args = parse_command_line_arguments()

    # Set up logging:
    initialize_logger(__name__, args.log_file)
    logger = logging.getLogger(__name__)

    main(
        args.output_file,
        args.ontoma_cache_dir,
    )
