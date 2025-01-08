"""Parser to build disease/target evidence based on Uniprot literature data."""

from __future__ import annotations

import argparse
import logging

from pyspark.sql import SparkSession
from pyspark.sql import functions as f
from uniprot_shared import UniprotShared

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)


class UniprotLiteratureExtractor(UniprotShared):
    """Class to extract Uniprot literature annotation from Unirot SPARQL API."""

    # SPARQL query to extract Uniprot literature data:
    UNIPROT_SPARQL_QUERY = """
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT
            ?protein 
            ?diseaseCrossrefs
            ?diseaseLabel
            ?diseaseComment
            (GROUP_CONCAT(DISTINCT ?source; SEPARATOR=", ") AS ?sources)
        
        WHERE
        {
            # Defining the universe of proteins:
            # BIND(<http://purl.uniprot.org/uniprot/B1AK53> AS ?protein) .
            ?protein a up:Protein ;
                up:organism taxon:9606 .

            ?protein rdfs:seeAlso ?diseaseCrossrefs .
            ?diseaseCrossrefs up:database <http://purl.uniprot.org/database/MIM> .
            ?diseaseCrossrefs rdfs:comment ?text .
            FILTER(CONTAINS(str(?text), "phenotype")) .
            ?disease_db rdfs:seeAlso ?diseaseCrossrefs .
            ?disease_db rdf:type up:Disease .
            ?disease_db skos:prefLabel ?diseaseLabel .
            ?disease_db rdfs:comment ?diseaseComment .
            OPTIONAL {
                ?linkToEvidence rdf:object ?disease_db ;
                                up:attribution ?attribution .
                ?attribution up:source ?source .
                ?source a up:Journal_Citation .
            }

            # Extracting gene to disease confidence:
            ?subject up:disease ?related_disease .
            ?subject rdfs:comment ?geneToDiseaseComment .

        }

        GROUP BY 
            ?protein 
            ?diseaseCrossrefs
            ?diseaseLabel
            ?diseaseComment
            ?geneToDiseaseComment
    """

    # List of columns for Uniprot literature evidence:
    EVIDENCE_COLUMNS = [
        "confidence",
        "datasourceId",
        "datatypeId",
        "diseaseFromSource",
        "diseaseFromSourceId",
        "diseaseFromSourceMappedId",
        "literature",
        "targetFromSourceId",
        "targetModulation",
    ]

    def __init__(self: UniprotLiteratureExtractor, spark: SparkSession) -> None:
        """Initialize the UniprotLiteratureExtractor class.

        Args:
            spark (SparkSession): The Spark session.
        """
        self.SPARK_SESSION = spark

    def extract_evidence_from_uniprot(
        self: UniprotLiteratureExtractor,
    ) -> UniprotLiteratureExtractor:
        """Extract Uniprot literature data from the Uniprot SPARQL API.

        Returns:
            UniprotLiteratureExtractor: The Uniprot literature extractor instance.
        """
        # Get data from Uniprot:
        literature_data = self.extract_uniprot_data()

        # Convert to a Spark DataFrame:
        self.evidence_dataframe = self.SPARK_SESSION.createDataFrame(
            literature_data
        ).select(
            # Extract disease information:
            f.col("diseaseLabel").alias("diseaseFromSource"),
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
            # Mapping geneToDiseaseComment to confidence:
            self.map_confidence(f.col("geneToDiseaseComment")).alias("confidence"),
            f.col("diseaseComment").alias("diseaseComment"),
            # Add default values:
            f.lit("uniprot_literature").alias("datasourceId"),
            f.lit("genetic_literature").alias("datatypeId"),
            f.lit("up_or_down").alias("targetModulation"),
        )

        return self


def main(
    output_file: str,
    ontoma_cache_dir: str,
) -> None:
    """Main function to extract Uniprot evidence.

    Args:
        output_file (str): The path to the output file.
        ontoma_cache_dir (str): The path to the OnToma cache directory.
    """
    # Get logger:
    logger.info("Starting Uniprot evidence parser.")
    logger.info(f"Output file: {output_file}")

    # Initialise spark session:
    spark = initialize_sparksession()

    # Extracting Uniprot evidence:
    uniprot_literature_evidence = (
        # Initialising the UniprotVariantsParser:
        UniprotLiteratureExtractor(spark)
        # Extract raw variant data:
        .extract_evidence_from_uniprot()
        # Map EFO terms:
        .add_efo_mapping(ontoma_cache_dir)
        # Accessing evidence data:
        .get_evidence(debug=True)
    )

    # Write data:
    logger.info("Writing data.")

    write_evidence_strings(
        uniprot_literature_evidence,
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
