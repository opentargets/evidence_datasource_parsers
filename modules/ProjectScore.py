"""Disease to target evidence parser for Project Score v2."""
import argparse
import logging

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as f
from pyspark.sql import types as t

from common.evidence import (
    GenerateDiseaseCellLines,
    initialize_logger,
    initialize_sparksession,
    write_evidence_strings,
)

# This parser is specific for the second release of Porject Score, so the publication identifier is hardcoded:
PMID = "38215750"


def generate_project_score_evidence(
    project_score_cell_lines: DataFrame, project_score_hits: DataFrame
) -> DataFrame:
    """Generate evidence strings for Project Score v2.

    Args:
        project_score_cell_lines (DataFrame): DataFrame containing cell line annotations
        project_score_hits (DataFrame): DataFrame containing evidence data

    Returns:
        DataFrame: DataFrame containing evidence strings
    """
    return (
        project_score_hits.select(
            f.col("targetSymbol").alias("targetFromSource"),
            f.col("diseaseName").alias("diseaseFromSource"),
            f.col("diseaseId").alias("diseaseFromSourceMappedId"),
            f.col("PRIORITY").cast(t.FloatType()).alias("resourceScore"),
            f.col("targetId").alias("targetFromSourceId"),
            f.array(f.lit(PMID)).alias("literature"),
            f.lit("crispr").alias("datasourceId"),
            f.lit("affected_pathway").alias("datatypeId"),
            f.lower(f.col("cancerType")).alias("cancerType"),
        )
        .join(project_score_cell_lines, on="cancerType", how="left")
        # Dropping the pancancer genes:
        .filter(f.col("cancerType") != "pancancer")
        # Cleaning table:
        .drop("cancerType")
    )


def main(
    spark: SparkSession,
    evid_file: str,
    cell_line_file: str,
    cell_passport_file: str,
    cell_line_to_uberon_mapping: str,
    out_file: str,
) -> None:
    """Main function for parsing Project Score v2 data.

    Args:
        spark (SparkSession):
        evid_file (str): File containing project score hits.
        cell_line_file (str): File containing project score cell lines.
        cell_passport_file (str): File containing Sanger Cell Model Passport data.
        cell_line_to_uberon_mapping (str): File containing tissue label to UBERON mapping.
        out_file (str): Output file name (gizpped json).
    """
    # Get logger:
    logger = logging.getLogger(__name__)

    # Logging command parameters:
    logger.info("Running Project Score parser with parameters:")
    logger.info(f"Evidence file: {evid_file}")
    logger.info(f"Cell types file: {cell_line_file}")
    logger.info(f"Cell passport file: {cell_passport_file}")
    logger.info(f"Cell line to uberon mapping: {cell_line_to_uberon_mapping}")
    logger.info(f"Output file: {out_file}")

    # Extract disease cell line data from cell passport file:
    cell_passport_data = GenerateDiseaseCellLines(
        spark, cell_passport_file, cell_line_to_uberon_mapping
    )

    passport_disease_cell_lines = cell_passport_data.generate_disease_cell_lines()

    # Joining disease cell-lines dataframe for Project Score:
    disease_cell_lines = (
        spark.read.csv(cell_line_file, sep="\t", header=True)
        .select(
            f.lower(f.col("CANCER_TYPE")).alias("cancerType"),
            f.col("CMP_ID").alias("id"),
        )
        .join(passport_disease_cell_lines, on="id", how="right")
        .groupBy("cancerType")
        .agg(f.collect_set("diseaseCellLine").alias("diseaseCellLines"))
    )
    # Are there any cancer types that are not in the cell line file?
    missing_cancer_types = disease_cell_lines.filter(f.col("diseaseCellLines").isNull())
    if missing_cancer_types.count() > 0:
        logger.warning(
            f"The following cancer types are not in the cell line file: {missing_cancer_types.collect()}"
        )
    else:
        logger.info("All cancer types are in the cell line file.")

    # Read gene based data and generate evidence strings:
    evidence_table = spark.read.csv(evid_file, sep="\t", header=True)

    project_score_evidence = generate_project_score_evidence(
        disease_cell_lines, evidence_table
    )
    logger.info(f"Number of evidence strings: {project_score_evidence.count()}")
    logger.info(
        f'Number of targets in evidence data: {project_score_evidence.select("targetFromSourceId").distinct().count()}'
    )
    logger.info(
        f'Number of diseases in evidence data: {project_score_evidence.select("diseaseFromSourceMappedId").distinct().count()}'
    )
    logger.info("Writing evidence strings to file...")

    write_evidence_strings(project_score_evidence, out_file)


def parse_args() -> argparse.Namespace:
    """Parsing command line arguments.

    Returns:
        argparse.Namespace: _description_
    """
    parser = argparse.ArgumentParser(
        description="Parse datasource parsers from Project Score"
    )

    parser.add_argument(
        "--evidence_file",
        help="Name of tsv file with the priority scores.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cell_types_file",
        help="Name of tsv file with applied cell line names and model passport identifier.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cell_passport_file",
        help="Csv file with the Sanger Cell Passport model data.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cell_line_to_uberon_mapping",
        help="Path to the manual curation file for mapping cell lines to uberon terms.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_file",
        help="Name of evidence file (gzip compressed json)",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--log_file",
        help="Name of log file.",
        type=str,
        required=False,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Parse arguments:
    args = parse_args()

    # Initialize logger:
    initialize_logger(__name__, args.log_file)

    # Initialize spark session:
    spark = initialize_sparksession()

    main(
        spark,
        args.evidence_file,
        args.cell_types_file,
        args.cell_passport_file,
        args.cell_line_to_uberon_mapping,
        args.output_file,
    )
