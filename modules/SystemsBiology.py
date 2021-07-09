import math
import argparse
import logging
import sys

import pandas as pd

def renormalize(n, start_range, new_range=[0.5, 1]):
    """
    A function to scale a value from a given range to a new range.

    Apply the function f(x) to n using and old (start_range) and a new range
    where f(x) = (dNewRange / dOldRange * (n - old_range_lower_bound)) + new_lower
    """

    delta1 = start_range[1] - start_range[0]
    delta2 = new_range[1] - new_range[0]

    max_new_range = max(new_range)
    min_new_range = min(new_range)

    if delta1 or delta2:
        try:
            normalized = (delta2 * (n - start_range[0]) / delta1) + new_range[0]
        except ZeroDivisionError:
            normalized = new_range[0]
    else:
        normalized = n

    # The formula results in values slightly smaller and larger than the boundaries of the new range
    if normalized > max_new_range:
        return max_new_range

    elif normalized < min_new_range:
        return min_new_range

    return round(normalized, 4)


def generate_score(row):
    """
    Score generation depends on the score type.
    """
    score_type = row['score_type']
    score = row['score']
    min_score = row['min_score']
    max_score = row['max_score']

    if score_type == 'p-value':
        parsed_score = renormalize(math.log10(score), [math.log10(max_score), math.log10(min_score)])
    elif score_type == 'rank':
        parsed_score = renormalize(score, [min_score, max_score])
    else:
        parsed_score = 0.75

    return parsed_score


def main(studyFile, evidenceFile, out_file):

    logging.info(f'Output file: {out_file}')

    # Reading evidence:
    logging.info(f'Evidence file: {evidenceFile}')
    evidence_df = pd.read_csv(evidenceFile, sep='\t')
    logging.info(f'Number of evidence: {len(evidence_df)}')
    logging.info(f'Number of target: {len(evidence_df.target_id.unique())}')
    logging.info(f'Number of disease: {len(evidence_df.disease_id.unique())}')

    # Reading study file:
    logging.info(f'Study description file: {studyFile}')
    publication_df = pd.read_csv(studyFile, sep='\t')
    logging.info(f'Number of studies: {len(publication_df)}')

    # Merging publication with evidence data:
    merged = evidence_df.merge(publication_df.drop('pmid', axis=1), on='gene_set_name', how='outer', indicator=True)

    # Checking if merging worked just fine:
    if len(merged.loc[merged._merge != 'both']) != 0:
        logging.warning(f'{len(merged.loc[merged._merge != "both"])} rows could not be joined.')
        logging.warning(merged.loc[merged._merge != "both"])

    # Generate evidence:
    merged = (
        merged
        .assign(
            diseaseFromSourceMappedId=merged.disease_id.apply(lambda x: x.split('/')[-1]),
            datasourceId='sysbio',
            datatypeId='affected_pathway',
            literature=merged.pmid.apply(lambda x: [str(x)]),
            pathways=merged.gene_set_name.apply(lambda x: [{'name': x}]),
            resourceScore=merged.apply(generate_score, axis=1)
        )
        .rename(columns={
            'target_id': 'targetFromSourceId',
            'disease_name': 'diseaseFromSource',
            'method': 'studyOverview'
        })
        .drop(['_merge', 'max_score', 'min_score', 'score_type', 'score', 'disease_id', 'pmid', 'gene_set_name'], axis=1)
        .to_json(out_file, compression='gzip', orient='records', lines=True)
    )

    logging.info('Evidence generation finished.')


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Generating target/disease evidence strings from SystemsBiology data.')
    parser.add_argument('-e', '--evidenceFile',
                        help='Name of tsv file with the per gene evidence file.',
                        type=str, required=True)
    parser.add_argument('-s', '--studyFile',
                        help='Name of tsv file with the study descriptions.',
                        type=str, required=True)
    parser.add_argument('-o', '--outputFile',
                        help='Name of gzip compressed json output file.',
                        type=str, required=True)
    parser.add_argument('-l', '--logFile',
                        help='Name of log file. If not provided logs are written to standard error.',
                        type=str, required=False)

    args = parser.parse_args()
    studyFile = args.studyFile
    evidenceFile = args.evidenceFile
    out_file = args.outputFile

    # If no logfile is specified, logs are written to stderr
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    main(studyFile, evidenceFile, out_file)
