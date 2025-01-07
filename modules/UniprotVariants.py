"""Parser to build disease/target evidence based on Uniprot variation data."""

from __future__ import annotations

import argparse
import logging

from pyspark.sql import DataFrame, SparkSession, Window
from pyspark.sql import functions as f
from uniprot_shared import UniprotShared

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping
from common.variant_rsid_mapping import RsIdMapper


class UniprotVariantsParser(UniprotShared):
    """Class to extract Uniprot variants from the Uniprot SPARQL API."""

    UNIPROT_SPARQL_QUERY = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT 
            ?protein 
            ?annotation 
            ?comment 
            ?variantRsId 
            ?diseaseLabel 
            ?diseaseComment 
            ?diseaseCrossrefs
            # Grouping pubmed links to a single string:
            (GROUP_CONCAT(DISTINCT ?source; SEPARATOR=", ") AS ?sources)
            ?referenceSequence
            ?begin
            ?end
            ?substitution
            ?geneToDiseaseComment
        WHERE
        {
            # Defining the universe of proteins:
            # BIND(<http://purl.uniprot.org/uniprot/B1AK53> AS ?protein) .
            ?protein a up:Protein ;
               up:organism taxon:9606 .

            ?protein up:annotation ?annotation .

            # Filtering the annotations to only include natural variants:
            ?annotation a up:Natural_Variant_Annotation . 

            ?annotation rdfs:comment ?comment .

            # Extracting the rsid of the variant if present:
            OPTIONAL { ?annotation rdfs:seeAlso ?variantRsId .}

            # Extracting substituted amino acid if present:
            OPTIONAL{ ?annotation up:substitution ?substitution .}

            # Extracting the disease annotations related to the natural variant, expected. Not interested in variants without disease annotations:
            ?annotation skos:related ?related_disease .
            ?related_disease skos:prefLabel ?diseaseLabel .
            ?related_disease rdfs:comment ?diseaseComment .
            ?related_disease rdfs:seeAlso ?diseaseCrossrefs .
            
            # Extracting gene to disease confidence:
            ?subject up:disease ?related_disease .
            ?subject rdfs:comment ?geneToDiseaseComment .

            # Extracting disease source if source is from OMIM:
            # FILTER(CONTAINS(STR(?diseaseCrossrefs), "mim"))
            ?diseaseCrossrefs up:database <http://purl.uniprot.org/database/MIM>

            # Extracting evidence links to the variant/disease relationship:
            OPTIONAL {
                ?linkToEvidence rdf:object ?annotation ;
                    up:attribution ?attribution .
                ?attribution up:source ?source .
                ?source a up:Journal_Citation .
            }

            # Extracting the sequence range of the variant (needs to be optional, as not all variants have a range): 
            OPTIONAL {
                ?annotation up:range ?range .
                ?range faldo:begin ?beginNode .
                ?range faldo:end ?endNode.
                ?endNode faldo:position ?end .
                ?beginNode faldo:position ?begin .
                ?beginNode faldo:reference ?sequence .
                ?sequence rdf:value ?value .
                BIND (substr(?value, ?begin, 1) as ?referenceSequence) .
            }
        }
        GROUP BY 
            ?protein 
            ?annotation 
            ?comment 
            ?variantRsId 
            ?substitution 
            ?diseaseLabel 
            ?diseaseComment 
            ?diseaseCrossrefs
            ?sequence
            ?begin
            ?end
            ?value
            ?referenceSequence
            ?geneToDiseaseComment
    """

    def __init__(
        self: UniprotVariantsParser,
        spark: SparkSession,
        ontoma_cache_dir: str,
    ):
        """Initialise the UniprotVariantsParser.

        Args:
            spark (SparkSession): The Spark session.
            rsid_cache (str): The rsid cache file.
            ontoma_cache_dir (str): The OnToma cache directory.
        """
        super().__init__(spark)
        self.ontoma_cache_dir = ontoma_cache_dir

    def extract_evidence_from_uniprot(
        self: UniprotVariantsParser,
    ) -> UniprotVariantsParser:
        """UniprotVariant specific pre-processing of the raw data."""
        variant_data = self.extract_uniprot_data()

        self.evidence_dataframe = self.SPARK_SESSION.createDataFrame(
            variant_data
        ).select(
            # Extract disease information:
            f.col("diseaseLabel").alias("diseaseFromSource"),
            f.concat(
                f.lit("OMIM:"),
                f.split("diseaseCrossrefs", "/").getItem(
                    f.size(f.split("diseaseCrossrefs", "/")) - 1
                ),
            ).alias("diseaseFromSourceId"),
            f.col("diseaseComment").alias("targetToDiseaseAnnotation"),
            # Extract target annotation
            f.col("protein").alias("targetFromSourceId"),
            f.col("comment").alias("variantToTargetAnnotation"),
            # Extract variant information - might be empty string:
            f.when(f.col("variantRsId") != "", f.col("variantRsId")).alias(
                "variantRsId"
            ),
            f.when(f.col("substitution") != "", f.col("substitution")).alias(
                "altAminoAcid"
            ),
            f.col("begin").alias("proteinPosition"),
            f.col("referenceSequence").alias("refAminoAcid"),
            # Extract pmids from sources:
            f.transform(
                f.split(f.col("sources"), ", "),
                lambda uri: f.split(uri, "/").getItem(f.size(f.split(uri, "/")) - 1),
            ).alias("literature"),
            "geneToDiseaseComment",
            # Mapping geneToDiseaseComment to confidence:
            self.map_confidence(f.col("geneToDiseaseComment")).alias("confidence"),
        )

        return self

    def map_rsids(
        self: UniprotVariantsParser, rsid_mapper: RsIdMapper
    ) -> UniprotVariantsParser:
        """Mapping rsids to variant ids.

        Args:
            rsid_mapper (RsIdMapper): rsid mapper object.

        Returns:
            UniprotVariantsParser: The UniprotVariantsParser.
        """
        # Get a list of unique rsids found in the dataset:
        unique_rsids = [
            row["variantRsId"]
            for row in self.evidence_dataframe.select("variantRsId")
            .distinct()
            .collect()
        ]

        logger.info(f"Number of unique rsids: {len(unique_rsids)}.")

        # Mapping rsids to variant ids:
        rsid_mapper.map_rsids(unique_rsids)

        # Save the cache file:
        rsid_mapper.update_cache_file()

        # Resolve variant ids with uniprot evidence:
        logger.info("Resolving variant ids.")
        self.evidence_dataframe = self.resolve_variant_ids(
            self.evidence_dataframe, rsid_mapper.get_mapped_variants()
        )

        return self

    @staticmethod
    def resolve_variant_ids(
        uniprot_variants: DataFrame, mapped_variants: DataFrame
    ) -> DataFrame:
        """Resolving variant ids for the Uniprot evidence.

        This function tries to find the best matching variantId for the
        Uniprot evidence based on the variantRsId and the predicted amino acid canges.

        Args:
            uniprot_variants (DataFrame): The Uniprot variants DataFrame.
            mapped_variants (DataFrame): The mapped variants DataFrame.

        Returns:
            DataFrame: The resolved Uniprot evidence DataFrame.
        """
        # Variant mappings are aggregated by rsId, disease and target to allow for selection of the correct variantId:
        window = Window.partitionBy(
            "diseaseFromSource", "targetFromSourceId", "variantRsId"
        )
        uniprot_columns = uniprot_variants.columns
        return (
            uniprot_variants
            # Join the uniprot variants with the variant mappings:
            .join(mapped_variants, on=["variantRsId", "targetFromSourceId"], how="left")
            # Create a boolean indicating if the variantId is accepted for the evidence:
            # This logic can be further refined to include more complex logic
            .withColumn(
                "isMatch",
                # Automatically accept biallelic variants:
                f.when(
                    ~f.col("isMultiAllelic"),
                    f.lit(True),
                    # Only accept multiallelic variants if they match the alternative amino acid matches the VEP annotation:
                )
                .when(f.col("altAminoAcid") == f.col("vepAltAminoAcid"), f.lit(True))
                .otherwise(f.lit(False)),
            )
            .withColumn(
                "isFailedWindow",
                # Flagging evidence/variantId windows where no match was found:
                f.when(
                    ~f.array_contains(f.collect_set("isMatch").over(window), True),
                    f.lit(True),
                ).otherwise(f.lit(False)),
            )
            # Accepting only matching variantids or failed windows:
            .filter(f.col("isFailedWindow") | f.col("isMatch"))
            .select(
                *uniprot_columns,
                f.col("variantId"),
                f.col("vepAltAminoAcid"),
                f.col("isMultiAllelic"),
            )
            .distinct()
        )

    def add_efo_mapping(
        self: UniprotVariantsParser, ontoma_cache_dir: str
    ) -> UniprotVariantsParser:
        """Add EFO mappings to the Uniprot evidence.

        Args:
            ontoma_cache_dir (str): The OnToma cache directory.

        Returns:
            UniprotVariantsParser: The UniprotVariantsParser.
        """
        self.evidence_dataframe = add_efo_mapping(
            self.evidence_dataframe, self.SPARK_SESSION, ontoma_cache_dir
        )
        # Add EFO mappings:
        logger.info("Adding EFO mappings.")
        self.evidence_dataframe = add_efo_mapping(
            self.evidence_dataframe, self.SPARK_SESSION, ontoma_cache_dir
        )

        return self

    def get_evidence(self: UniprotVariantsParser, debug: bool = False) -> DataFrame:
        """Get the Uniprot evidence.

        Returns:
            DataFrame: The Uniprot evidence.
        """
        # Return all columns for debugging purposes:
        if debug:
            return self.evidence_dataframe

        return self.evidence_dataframe.select(
            f.lit("uniprot_variants").alias("datasourceId"),
            f.lit("genetic_association").alias("datatypeId"),
            f.lit("up_or_down").alias("targetModulation"),
        )


def main(
    rsid_cache: str,
    output_file: str,
    ontoma_cache_dir: str,
) -> None:
    # Get logger:
    logger.info("Starting Uniprot evidence parser.")
    logger.info(f"Output file: {output_file}")

    # Initialise spark session:
    spark = initialize_sparksession()

    # Initialising rsid mapper with the provided cache file:
    rsid_mapper = RsIdMapper(spark, rsid_cache)

    # Extracting Uniprot evidence:
    uniprot_variants_evidence = (
        # Initialising the UniprotVariantsParser:
        UniprotVariantsParser(spark, ontoma_cache_dir)
        # Extract raw variant data:
        .extract_evidence_from_uniprot()
        # map rsIDs to variant IDs:
        .map_rsids(rsid_mapper)
        # Map EFO terms:
        .add_efo_mapping(ontoma_cache_dir)
        # Accessing evidence data:
        .get_evidence(debug=True)
    )

    # Write data:
    logger.info("Writing data.")

    write_evidence_strings(
        uniprot_variants_evidence,
        output_file,
    )


def parse_command_line_arguments() -> argparse.Namespace:
    """Parsing command line arguments for the Uniprot evidence parser.

    Returns:
        argparse.Namespace: The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Parser to build disease/target evidence based on Uniprot variation data."
    )
    parser.add_argument(
        "--rsid_cache",
        type=str,
        help="RsId cache file containing rsid to variant id mapping.",
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
        args.rsid_cache,
        args.output_file,
        args.ontoma_cache_dir,
    )
