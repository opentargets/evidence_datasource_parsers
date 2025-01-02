"""Parser to build disease/target evidence based on Uniprot variation data."""

from __future__ import annotations

import argparse
import logging
import tempfile
from typing import Any

import requests
from pyspark.sql import DataFrame, SparkSession, Window
from pyspark.sql import functions as f
from pyspark.sql import types as t
from pyspark.sql.utils import AnalysisException
from SPARQLWrapper import JSON, SPARQLWrapper

from common.evidence import (
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)
from common.ontology import add_efo_mapping


class UniprotVariantsExtractor:
    """Class to extract Uniprot variants from the Uniprot SPARQL API."""

    UNIPROT_SPARQL_ENDPOINT = "https://sparql.uniprot.org/sparql"

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

    @classmethod
    def extract_variants_data(
        cls: type[UniprotVariantsExtractor],
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

        variant_data = []
        for result in results["results"]["bindings"]:
            variant = {}
            for key, value in result.items():
                # Extracting only MIM disease cross-references:
                if key == "disease_crossrefs" and "mim" not in value["value"]:
                    continue

                if value["type"] == "literal":
                    variant[key] = value["value"]
                elif "mim" in value["value"]:
                    variant[key] = value["value"]
                elif value["type"] == "uri":
                    variant[key] = cls.get_uri_leaf(value["value"])

            variant_data.append(variant)

        logger.info(f"Number of disease/target/variant evidence: {len(variant_data)}.")
        return variant_data

    @classmethod
    def get_dataframe(
        cls: type[UniprotVariantsExtractor], spark: SparkSession
    ) -> DataFrame:
        """Convert VEP data to a Spark DataFrame and minimally process it.

        Args:
            spark (SparkSession): The Spark session.

        Returns:
            DataFrame: The processed VEP data.
        """
        variant_data = cls.extract_variants_data()

        return spark.createDataFrame(variant_data).select(
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
        )

    @staticmethod
    def get_uri_leaf(uri: str) -> str:
        return uri.split("/")[-1]


class rsidMapper:
    CACHE_SCHEMA = t.StructType(
        [
            t.StructField("variantRsId", t.StringType(), True),
            t.StructField("variantId", t.StringType(), True),
            t.StructField("vepProteinPosition", t.IntegerType(), True),
            t.StructField("vepRefAminoAcid", t.StringType(), True),
            t.StructField("vepAltAminoAcid", t.StringType(), True),
            t.StructField("targetFromSourceId", t.StringType(), True),
            t.StructField("vcf_string", t.ArrayType(t.StringType(), True), True),
            t.StructField("isMultiAllelic", t.BooleanType(), True),
        ]
    )

    API_RESPONSE_SCHEMA = t.StructType(
        [
            t.StructField("start", t.LongType(), nullable=False),
            t.StructField("end", t.LongType(), nullable=False),
            t.StructField("input", t.StringType(), nullable=False),
            t.StructField("allele_string", t.StringType(), nullable=False),
            t.StructField(
                "transcript_consequences",
                t.ArrayType(
                    t.StructType(
                        [
                            t.StructField(
                                "protein_start", t.IntegerType(), nullable=True
                            ),
                            t.StructField(
                                "variant_allele", t.StringType(), nullable=False
                            ),
                            t.StructField(
                                "swissprot", t.ArrayType(t.StringType()), nullable=True
                            ),
                            t.StructField("amino_acids", t.StringType(), nullable=True),
                        ]
                    )
                ),
                nullable=True,
            ),
            t.StructField("vcf_string", t.ArrayType(t.StringType()), nullable=True),
            t.StructField("seq_region_name", t.StringType(), nullable=False),
        ]
    )

    API_SIZE_LIMIT = 200

    API_URL = "https://rest.ensembl.org/vep/human/id"
    REQUEST_HEADERS = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    RESQUEST_PARAMS = {
        "uniprot": 1,
        "vcf_string": 1,
    }

    def __init__(self, spark: SparkSession, rsid_cache: str) -> None:
        # Store the cache file name:
        self.rsid_cache = rsid_cache
        self.spark = spark

        # Try to load the cache file:
        try:
            logger.info(f"Reading cache from: {rsid_cache}.")
            self.rsid_cache_df = spark.read.parquet(rsid_cache)
        except AnalysisException:
            logger.info(
                f"The provided cache file could not be read. Creating new cache here: {rsid_cache}."
            )
            self.rsid_cache_df = spark.createDataFrame([], schema=self.CACHE_SCHEMA)

        # Get a list of already mapped rsids:
        self.cached_rsids: list[str] = [
            row["variantRsId"]
            for row in self.rsid_cache_df.select("variantRsId").distinct().collect()
        ]
        # Register with a temporary view:
        self.rsid_cache_df.createOrReplaceTempView("cached_rsids")

        logger.info(f"Number of cached rsids: {len(self.cached_rsids)}.")

    def update_cache_file(self) -> None:
        """Save the cache file in the provided location."""
        # Create a temporary file
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".parquet")
        temp_file_path = temp_file.name
        # Write the updated DataFrame to the temporary Parquet file
        self.rsid_cache_df.write.mode("overwrite").parquet(temp_file_path)

        (
            # Read the temporary Parquet file
            self.spark.read.parquet(temp_file_path)
            .write.mode("overwrite")
            .parquet(self.rsid_cache)
        )
        temp_file.close()
        # Refresh the table to invalidate the cache
        self.spark.sql("REFRESH TABLE cached_rsids")

    def get_mapped_variants(self) -> DataFrame:
        """Return the cache DataFrame.

        Returns:
            DataFrame: The cache DataFrame.
        """
        return self.rsid_cache_df

    def map_rsids(self, rsids) -> None:
        """Map rsids to variant ids. Adds the mapping to the cache.

        Args:
            rsids (list[str]): List of rsids to map.
        """
        # Get the missing rsids:
        missing_rsids = [rsid for rsid in rsids if rsid not in self.cached_rsids]

        # Chunking the missing rsids into accepted chunks:
        rsid_chunks = [
            missing_rsids[i : i + self.API_SIZE_LIMIT]
            for i in range(0, len(missing_rsids), self.API_SIZE_LIMIT)
        ]

        # Report missingness:
        logger.info(
            f"Number of missing rsids that needs to be mapped: {len(missing_rsids)}."
        )
        logger.info(f"Missing rsids are broken down into {len(rsid_chunks)} chunks.")

        # Mapping the missing rsids:
        for index, rsid_chunk in enumerate(rsid_chunks):
            logger.info(f"Mapping chunk {index}...")
            # Get mapping for the chunk:
            mapped_chunk = self._map_rsids_chunk(rsid_chunk)
            # Adding chunk to the cache:
            self.rsid_cache_df = self.rsid_cache_df.unionByName(mapped_chunk)

        logger.info("Mapping rsids completed.")

    def _map_rsids_chunk(self, rsids) -> DataFrame:
        vep_annotated_data = self._get_vep_for_rsid(rsids)

        # Ensure the vcf_string is a list - neccesary as the VEP schema is inconsistent:
        for data in vep_annotated_data:
            if not isinstance(data["vcf_string"], list):
                data["vcf_string"] = [data["vcf_string"]]

        return (
            # Read VEP output into a DataFrame following a simplified response schema:
            self.spark.createDataFrame(vep_annotated_data, self.API_RESPONSE_SCHEMA)
            # Explode transcript consequences and uniprot ids:
            .withColumn("transcript_consequences", f.explode("transcript_consequences"))
            .select("*", "transcript_consequences.*")
            .withColumn("uniprot_id", f.explode("swissprot"))
            # Do transformations to get the final schema:
            .select(
                # Selecting only the necessary columns:
                f.col("input").alias("variantRsId"),
                # TODO: The variant id should be generated from the vcf_string field:
                f.concat_ws(
                    "_",
                    "seq_region_name",
                    "start",
                    f.split(f.col("allele_string"), "/")[0],
                    "variant_allele",
                ).alias("variantId"),
                # Extracting the protein position and substitued amino acids:
                f.col("protein_start").alias("vepProteinPosition"),
                f.split(f.col("amino_acids"), "/")[0].alias("vepRefAminoAcid"),
                f.when(
                    f.split(f.col("amino_acids"), "/")[1] != "*",
                    f.split(f.col("amino_acids"), "/")[1],
                ).alias("vepAltAminoAcid"),
                # Stripping uniprot version from id:
                f.split(f.col("uniprot_id"), r"\.")[0].alias("targetFromSourceId"),
                f.col("vcf_string"),
                # Flag multiallelic variants:
                (f.size("vcf_string") > 1).alias("isMultiAllelic"),
            )
            # Dropping non-protein coding consequeces:
            .filter(f.col("targetFromSourceId").isNotNull())
            .distinct()
            .orderBy("variantRsId", "variantId")
        )

    def _get_vep_for_rsid(self, variants_to_map) -> list:
        """Get VEP data for a list of variants.

        Args:
            variants_to_map (list[str]): List of variants to map.

        Returns:
            dict[str, Any]: The VEP data for the given variants.
        """
        response = requests.post(
            self.API_URL,
            headers=self.REQUEST_HEADERS,
            params=self.RESQUEST_PARAMS,
            json={"ids": variants_to_map},
        )

        return response.json()


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

    # Extracting Uniprot evidence:
    uniprot_variants = (
        UniprotVariantsExtractor.get_dataframe(spark)
        # For testing purposes, limit the number of rows:
        # .orderBy(f.rand())
        # .limit(20)
    )

    # Get unique rsids:
    unique_rsids = [
        row["variantRsId"]
        for row in uniprot_variants.select("variantRsId").distinct().collect()
    ]

    logger.info(f"Number of unique rsids: {len(unique_rsids)}.")

    # Initialising rsid mapper with the provided cache file:
    rsid_mapper = rsidMapper(spark, rsid_cache)

    # Mapping rsids to variant ids:
    rsid_mapper.map_rsids(unique_rsids)

    # Save the cache file:
    rsid_mapper.update_cache_file()

    # Resolve variant ids with uniprot evidence:
    logger.info("Resolving variant ids.")
    resolved_evidence = resolve_variant_ids(
        uniprot_variants, rsid_mapper.get_mapped_variants()
    )

    # Add EFO mappings:
    logger.info("Adding EFO mappings.")
    evidence_df = add_efo_mapping(resolved_evidence, spark, ontoma_cache_dir)

    # Write data:
    logger.info("Writing data.")
    evidence_df.write.mode("overwrite").parquet("uniprot_debug.parquet")
    # write_evidence_strings(
    #     (
    #         evidence_df
    #         # TODO: Based on targetToDiseaseAnnotation classify confidence
    #         # TODO: Based on variantToTargetAnnotation classify direction of effect/target modulation
    #         # The reference and altenate amino acids are there for troubleshooting purposes.
    #         .drop(
    #             "targetToDiseaseAnnotation",
    #             "variantToTargetAnnotation",
    #             "altAminoAcid",
    #             "proteinPosition",
    #             "refAminoAcid",
    #             "isMultiAllelic",
    #             "vepAltAminoAcid",
    #         ).withColumns(
    #             {
    #                 "datasourceId": f.lit("uniprot_variants"),
    #                 "datatypeId": f.lit("genetic_association"),
    #                 "confidence": f.lit("high"),
    #                 "targetModulation": f.lit("up_or_down"),
    #             }
    #         )
    #     ),
    #     output_file,
    # )


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
