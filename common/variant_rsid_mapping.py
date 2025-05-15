"""Tools to map variant rsids to variant location and variant identifiers.

Considerations:
- This tool uses the Ensembl Variant Effect Predictor (VEP) REST API to map rsids to variant locations.
- As the mapping can be time consuming, the tool uses a cache file to store the already mapped rsids.
- The cache file can be updated with new mappings once the mapping is complete.
"""

from __future__ import annotations

import logging
import tempfile

import requests
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t
from pyspark.sql.utils import AnalysisException

logger = logging.getLogger(__name__)


class RsIdMapper:
    """Class to map rsids to variant ids and locations."""

    # Schema for the cache file is required to initialise empty DataFrame:
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

    # Schema for the VEP API response:
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

    # The maximum size of rsIDs that can be submitted to the API at once:
    API_SIZE_LIMIT = 200

    API_URL = "https://rest.ensembl.org/vep/human/id"
    REQUEST_HEADERS = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    # Extra parameters for the API request:
    RESQUEST_PARAMS = {
        "uniprot": 1,
        "vcf_string": 1,
    }

    def __init__(self: RsIdMapper, spark: SparkSession, rsid_cache: str) -> None:
        """Initialize the rsidMapper class.

        As the mapping can be time consuming, the class uses a cache file to store the
        already mapped rsids. The cache file is used to avoid mapping the same rsids multiple times.
        If the cache file cannot be read, a new cache file is created and all rsids are mapped.

        Args:
            spark (SparkSession): The Spark session.
            rsid_cache (str): The location of the cache file.
        """
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

    def update_cache_file(self: RsIdMapper) -> None:
        """Save the cache file in the provided location.

        The cache file is saved in a temporary file and then moved to the provided location.

        Raises:
            AnalysisException: If the cache file cannot be read.
        """
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

    def map_rsids(self: RsIdMapper, rsids: list[str]) -> None:
        """Map rsids to variant ids. Adds the mapping to the cache.

        Args:
            rsids (list[str]): List of rsids to map.
        """
        # Get a unique list of rsIDs that are not already cached:
        missing_rsids = [
            rsid
            for rsid in set(rsids)
            if rsid not in self.cached_rsids and rsid is not None
        ]

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

    def _map_rsids_chunk(self: RsIdMapper, rsids: list[str]) -> DataFrame:
        """Map a chunk of rsids to variant ids.

        Args:
            rsids (list[str]): List of rsids to map.

        Returns:
            DataFrame: The mapped rsids as spark dataframe.
        """
        vep_annotated_data = self._get_vep_for_rsid(rsids)

        # Ensure the vcf_string is a list - neccesary as the VEP schema is inconsistent:
        for data in vep_annotated_data:
            if not isinstance(data, dict):
                logger.error(f"VEP response is not a dictionary! Found: {type(data)}.")
                raise ValueError(
                    f"VEP response is not a dictionary! Found: {type(data)}.\n{data}"
                )
            if not isinstance(data["vcf_string"], list):
                data["vcf_string"] = [data["vcf_string"]]

        return (
            # Read VEP output into a DataFrame following a simplified response schema:
            self.spark.createDataFrame(vep_annotated_data, self.API_RESPONSE_SCHEMA)
            # Explode transcript consequences and uniprot ids:
            .withColumn(
                "transcript_consequences", f.explode_outer("transcript_consequences")
            )
            .select("*", "transcript_consequences.*")
            .withColumn("uniprot_id", f.explode_outer("swissprot"))
            # Dropping mappings from non-canonical chromosomes:
            .filter(f.col("seq_region_name").rlike(r"^\d{1,2}|X|Y$"))
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
                (f.size(f.col("vcf_string")) > 1).alias("isMultiAllelic"),
            )
            # Dropping non-protein coding consequeces:
            # .filter(f.col("targetFromSourceId").isNotNull())
            .distinct()
            .orderBy("variantRsId", "variantId")
        )

    def _get_vep_for_rsid(self: RsIdMapper, variants_to_map: list) -> list:
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
        if response.status_code != 200:
            logger.warning(
                f"Failed to map rsids: {variants_to_map}. Response: {response.text}"
            )
            return []

        return response.json()
