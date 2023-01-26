#!/usr/bin/env python3
"""Evidence parser for the Genomics England PanelApp data."""

import argparse
import logging
import re

import requests
from pyspark.sql import functions as f, SparkSession

from common.ontology import add_efo_mapping
from common.evidence import write_evidence_strings


class PanelAppEvidenceGenerator:

    # Fixes which are applied to original PanelApp phenotype strings *before* splitting them by semicolon.
    PHENOTYPE_BEFORE_SPLIT_RE = {
        # Fixes for specific records.
        r"\(HP:0006574;\);": r"(HP:0006574);",
        r"Abruzzo-Erickson;syndrome": r"Abruzzo-Erickson syndrome",
        r"Deafness, autosomal recessive; 12": r"Deafness, autosomal recessive, 12",
        r"Waardenburg syndrome, type; 3": r"Waardenburg syndrome, type 3",
        r"Ectrodactyly, ectodermal dysplasia, and cleft lip/palate syndrome; 3": r"Ectrodactyly, ectodermal dysplasia, and cleft lip/palate syndrome, 3",
        # Remove curly braces. They are *sometimes* (not consistently) used to separate disease name and OMIM code, for
        # example: "{Bladder cancer, somatic}, 109800", and interfere with regular expressions for extraction.
        r"[{}]": r"",
        # Fix cases like "Aarskog-Scott syndrome, 305400Mental retardation, X-linked syndromic 16, 305400", where
        # several phenotypes are glued to each other due to formatting errors.
        r"(\d{6})([A-Za-z])": r"$1;$2",
        # Replace all tab/space sequences with a single space.
        r"[\t ]+": r" ",
        # Remove leading and/or trailing spaces around semicolons.
        r" ?; ?": r";",
    }

    # Cleanup regular expressions which are applied to the phenotypes *after* splitting.
    PHENOTYPE_AFTER_SPLIT_RE = (
        r" \(no OMIM number\)",
        r" \(NO phenotype number in OMIM\)",
        r"(no|No|NO) OMIM( phenotype|number|entry|NUMBER|NUMBER OR DISEASE)?",
        r"[( ]*(from )?PMID:? *\d+[ ).]*",
    )

    # Regular expressions for extracting ontology information.
    LEADING = r"[ ,-]*"
    SEPARATOR = r"[:_ #]*"
    TRAILING = r"[:.]*"
    OMIM_RE = LEADING + r"(OMIM|MIM)?" + SEPARATOR + r"(\d{6})" + TRAILING
    OTHER_RE = (
        LEADING
        + r"(OrphaNet: ORPHA|Orphanet|ORPHA|HP|MONDO)"
        + SEPARATOR
        + r"(\d+)"
        + TRAILING
    )

    # Regular expression for filtering out bad PMID records.
    PMID_FILTER_OUT_RE = r"224,614,752,030,146,000,000,000"

    # Regular expressions for extracting publication information from the API raw strings.
    PMID_RE = [
        (
            # Pattern 1, e.g. '15643612' or '28055140, 27333055, 23063529'.
            r"^"  # Start of the string.
            r"[\d, ]+"  # A sequence of digits, commas and spaces.
            r"(?: |$)"  # Ending either with a space, or with the end of the string.
        ),
        (
            # Pattern 2, e.g. '... observed in the patient. PMID: 1908107 - publication describing function of ...'
            r"(?:PubMed|PMID)"  # PubMed or a PMID prefix.
            r"[: ]*"  # An optional separator (spaces/colons).
            r"[\d, ]+"  # A sequence of digits, commas and spaces.
        ),
    ]

    def __init__(self, debug_output_prefix=None):
        self.spark = SparkSession.builder.appName("panelapp_parser").getOrCreate()
        self.debug_output_phenotypes_filename = None
        self.debug_output_pmids_filename = None
        self.debug_output_pmids_data = set()
        if debug_output_prefix:
            self.debug_output_phenotypes_filename = (
                f"{debug_output_prefix}.phenotypes.tsv"
            )
            self.debug_output_pmids_filename = f"{debug_output_prefix}.pmids.tsv"

    def __del__(self):
        if self.debug_output_pmids_filename:
            with open(self.debug_output_pmids_filename, "w") as outfile:
                outfile.write(
                    "Source PMID string\tExtracted PMIDs\tNumber of extracted PMIDs\n"
                )
                for line in sorted(set(self.debug_output_pmids_data)):
                    outfile.write(line)

    def generate_panelapp_evidence(
        self, input_file: str, output_file: str, cache_dir: str
    ) -> None:
        logging.info("Filter and extract the necessary columns.")
        panelapp_df = self.spark.read.csv(input_file, sep=r"\t", header=True)
        # Panel version can be either a single number (e.g. 1), or two numbers separated by a dot (e.g. 3.14). We cast
        # either representation to float to ensure correct filtering below. (Note that conversion to float would not
        # work in the general case, because 3.4 > 3.14, but we only need to compare relative to 1.0.)
        panelapp_df = panelapp_df.withColumn(
            "Panel Version",
            panelapp_df["Panel Version"].cast("float").alias("Panel Version"),
        )
        panelapp_df = (
            panelapp_df.filter(
                ((f.col("List") == "green") | (f.col("List") == "amber"))
                & (f.col("Panel Version") >= 1.0)
                & (f.col("Panel Status") == "PUBLIC")
            ).select(
                "Symbol",
                "Panel Id",
                "Panel Name",
                "List",
                "Mode of inheritance",
                "Phenotypes",
            )
            # The full original records are not redundant; however, uniqueness on a subset of fields is not guaranteed.
            .distinct()
        )

        logging.info(
            "Fix typos and formatting errors which would interfere with phenotype splitting."
        )
        panelapp_df = panelapp_df.withColumn("cleanedUpPhenotypes", f.col("Phenotypes"))
        for regexp, replacement in self.PHENOTYPE_BEFORE_SPLIT_RE.items():
            panelapp_df = panelapp_df.withColumn(
                "cleanedUpPhenotypes",
                f.regexp_replace(f.col("cleanedUpPhenotypes"), regexp, replacement),
            )

        logging.info("Split and explode the phenotypes.")
        panelapp_df = panelapp_df.withColumn(
            "cohortPhenotypes",
            f.array_distinct(f.split(f.col("cleanedUpPhenotypes"), ";")),
        ).withColumn("phenotype", f.explode(f.col("cohortPhenotypes")))

        logging.info(
            "Remove specific patterns and phrases which will interfere with ontology extraction and mapping."
        )
        panelapp_df = panelapp_df.withColumn("diseaseFromSource", f.col("phenotype"))
        for regexp in self.PHENOTYPE_AFTER_SPLIT_RE:
            panelapp_df = panelapp_df.withColumn(
                "diseaseFromSource",
                f.regexp_replace(f.col("diseaseFromSource"), f"({regexp})", ""),
            )

        logging.info(
            "Extract ontology information, clean up and filter the split phenotypes."
        )
        panelapp_df = (
            panelapp_df
            # Extract Orphanet/MONDO/HP ontology identifiers and remove them from the phenotype string.
            .withColumn(
                "ontology_namespace",
                f.regexp_extract(f.col("diseaseFromSource"), self.OTHER_RE, 1),
            )
            .withColumn(
                "ontology_namespace",
                f.regexp_replace(
                    f.col("ontology_namespace"), "OrphaNet: ORPHA", "ORPHA"
                ),
            )
            .withColumn(
                "ontology_id",
                f.regexp_extract(f.col("diseaseFromSource"), self.OTHER_RE, 2),
            )
            .withColumn(
                "ontology",
                f.when(
                    (f.col("ontology_namespace") != "") & (f.col("ontology_id") != ""),
                    f.concat(
                        f.col("ontology_namespace"), f.lit(":"), f.col("ontology_id")
                    ),
                ),
            )
            .withColumn(
                "diseaseFromSource",
                f.regexp_replace(f.col("diseaseFromSource"), f"({self.OTHER_RE})", ""),
            )
            # Extract OMIM identifiers and remove them from the phenotype string.
            .withColumn(
                "omim_id", f.regexp_extract(f.col("diseaseFromSource"), self.OMIM_RE, 2)
            )
            .withColumn(
                "omim",
                f.when(
                    f.col("omim_id") != "", f.concat(f.lit("OMIM:"), f.col("omim_id"))
                ),
            )
            .withColumn(
                "diseaseFromSource",
                f.regexp_replace(f.col("diseaseFromSource"), f"({self.OMIM_RE})", ""),
            )
            # Choose one of the ontology identifiers, keeping OMIM as a priority.
            .withColumn(
                "diseaseFromSourceId",
                f.when(f.col("omim").isNotNull(), f.col("omim")).otherwise(
                    f.col("ontology")
                ),
            )
            .drop("ontology_namespace", "ontology_id", "ontology", "omim_id", "omim")
            # Clean up the final split phenotypes.
            .withColumn(
                "diseaseFromSource",
                f.regexp_replace(f.col("diseaseFromSource"), r"\(\)", ""),
            )
            .withColumn("diseaseFromSource", f.trim(f.col("diseaseFromSource")))
            .withColumn(
                "diseaseFromSource",
                f.when(f.col("diseaseFromSource") != "", f.col("diseaseFromSource")),
            )
            # Remove low quality records, where the name of the phenotype string starts with a question mark.
            .filter(
                ~(
                    (f.col("diseaseFromSource").isNotNull())
                    & (f.col("diseaseFromSource").startswith("?"))
                )
            )
            # Remove duplication caused by cases where multiple phenotypes within the same record fail to generate any
            # phenotype string or ontology identifier.
            .distinct()
            # For records where we were unable to determine either a phenotype string nor an ontology identifier,
            # substitute the panel name instead.
            .withColumn(
                "diseaseFromSource",
                f.when(
                    (f.col("diseaseFromSource").isNull())
                    & (f.col("diseaseFromSourceId").isNull()),
                    f.col("Panel Name"),
                ).otherwise(f.col("diseaseFromSource")),
            )
            .persist()
        )

        logging.info("Fetch and join literature references.")
        all_panel_ids = panelapp_df.select("Panel Id").toPandas()["Panel Id"].unique()
        literature_references = self.fetch_literature_references(all_panel_ids)
        panelapp_df = panelapp_df.join(
            literature_references, on=["Panel Id", "Symbol"], how="left"
        )

        if self.debug_output_phenotypes_filename:
            logging.info("Output tables for debugging purposes, if requested.")
            (
                panelapp_df.select(
                    "Phenotypes",  # Original, unaltered string with all phenotypes.
                    "cleanedUpPhenotypes",  # String with phenotypes after pre-split cleanup.
                    "phenotype",  # Individual phenotype after splitting.
                    "diseaseFromSource",  # Final cleaned up disease name.
                    "diseaseFromSourceId",  # Final cleaned up disease ID.
                )
                .distinct()
                .toPandas()
                .to_csv(self.debug_output_phenotypes_filename, sep="\t", index=False)
            )

        logging.info(
            "Drop unnecessary fields and populate the final evidence string structure."
        )
        evidence_df = (
            panelapp_df.drop("Phenotypes", "cleanedUpPhenotypes", "phenotype")
            # allelicRequirements requires a list, but we always only have one value from PanelApp.
            .withColumn(
                "allelicRequirements",
                f.when(
                    f.col("Mode of inheritance").isNotNull(),
                    f.array(f.col("Mode of inheritance")),
                ),
            )
            .drop("Mode of inheritance")
            .withColumnRenamed("List", "confidence")
            .withColumn("datasourceId", f.lit("genomics_england"))
            .withColumn("datatypeId", f.lit("genetic_literature"))
            # diseaseFromSourceId populated above
            # literature populated above
            .withColumnRenamed("Panel Id", "studyId")
            .withColumnRenamed("Panel Name", "studyOverview")
            .withColumnRenamed("Symbol", "targetFromSourceId")
            # Some residual duplication is caused by slightly different representations from `cohortPhenotypes` being
            # cleaned up to the same representation in `diseaseFromSource`, for example "Pontocerebellar hypoplasia type
            # 2D (613811)" and "Pontocerebellar hypoplasia type 2D, 613811".
            .distinct()
        )

        evidence_df = add_efo_mapping(
            evidence_strings=evidence_df,
            spark_instance=self.spark,
            ontoma_cache_dir=cache_dir,
        )
        logging.info("Disease mappings have been added.")

        write_evidence_strings(evidence_df, output_file)
        logging.info(
            f"{evidence_df.count()} evidence strings have been saved to {output_file}"
        )

    def fetch_literature_references(self, all_panel_ids):
        """Queries the PanelApp API to extract all literature references for (panel ID, gene symbol) combinations."""
        publications = []  # Contains tuples of (panel ID, gene symbol, PubMed ID).
        for panel_id in all_panel_ids:
            logging.debug(f"Fetching literature references for panel {panel_id!r}.")
            url = f"https://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}"
            panel_data = requests.get(url).json()

            # The source data and the online data might not be in-sync, due to requesting retired panels.
            # We still keep these entries, but won't try to fetch gene data:
            if "gene" not in panel_data:
                logging.warn(f"Panel info could not retrieved for panel id: {panel_id}")
                continue

            for gene in panel_data["genes"]:
                for publication_string in gene["publications"]:
                    publications.extend(
                        [
                            (panel_id, gene["gene_data"]["gene_symbol"], pubmed_id)
                            for pubmed_id in self.extract_pubmed_ids(publication_string)
                        ]
                    )
        # Group by (panel ID, gene symbol) pairs and convert into a PySpark dataframe.
        return (
            self.spark.createDataFrame(
                publications, schema=("Panel ID", "Symbol", "literature")
            )
            .groupby(["Panel ID", "Symbol"])
            .agg(collect_set("literature").alias("literature"))
        )

    def extract_pubmed_ids(self, publication_string):
        """Parses the publication information from the PanelApp API and extracts PubMed IDs."""
        publication_string = (
            publication_string.encode("ascii", "ignore")
            .decode()  # To get rid of zero width spaces and other special characters.
            .strip()
            .replace("\n", "")
            .replace("\r", "")
        )
        pubmed_ids = []

        if not re.match(self.PMID_FILTER_OUT_RE, publication_string):
            for regexp in self.PMID_RE:  # For every known representation pattern...
                for occurrence in re.findall(
                    regexp, publication_string
                ):  # For every occurrence of this pattern...
                    pubmed_ids.extend(
                        re.findall(r"(\d+)", occurrence)
                    )  # Extract all digit sequences (PubMed IDs).

        # Filter out:
        # * 0 as a value, because it is a placeholder for a missing ID;
        # * PubMed IDs which are too short or too long.
        pubmed_ids = {
            pubmed_id
            for pubmed_id in pubmed_ids
            if pubmed_id != "0" and len(pubmed_id) <= 8
        }

        if self.debug_output_pmids_filename:
            # Output the original PMID string, the extracted PMIDs, and the total number of the records extracted.
            self.debug_output_pmids_data.add(
                f'{publication_string}\t{" ".join(sorted(pubmed_ids))}\t{len(pubmed_ids)}\n'
            )

        return pubmed_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-file",
        required=True,
        type=str,
        help="Input TSV file with the table containing association details.",
    )
    parser.add_argument(
        "--output-file",
        required=True,
        type=str,
        help="Output JSON file (gzip-compressed) with the evidence strings.",
    )
    parser.add_argument(
        "--debug-output-prefix",
        required=False,
        type=str,
        help="If specified, a number of files will be created with this prefix to store phenotype/PMID cleanup data "
        "for debugging purposes.",
    )
    parser.add_argument(
        "--cache_dir",
        required=False,
        help="Directory to store the OnToma cache files in.",
    )

    args = parser.parse_args()
    PanelAppEvidenceGenerator(args.debug_output_prefix).generate_panelapp_evidence(
        input_file=args.input_file,
        output_file=args.output_file,
        cache_dir=args.cache_dir,
    )
