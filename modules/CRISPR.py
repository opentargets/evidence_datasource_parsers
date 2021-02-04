import sys
import logging
import datetime
import pandas as pd
import argparse


# A few genes do not have Ensembl IDs in the data file provided
CRISPR_SYMBOL_MAPPING={
    'CASC5': 'ENSG00000137812',
    'CIRH1A': 'ENSG00000141076',
    'EFTUD1': 'ENSG00000140598',
    'ENSG00000163660': 'ENSG00000163660',
    'KIAA0947': 'ENSG00000164151',
    'KIAA1432': 'ENSG00000107036',
    'NDNL2': 'ENSG00000185115',
    'SRPR': 'ENSG00000182934',
    'ZNF259': 'ENSG00000109917'
}



def main():

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse essential cancer genes identified using CRISPR assays and prioritised in Project Score')
    parser.add_argument('-d', '--descriptions_file',
                        help='Name of tsv file with the description of the method per cancer type',
                        type=str, required=True)
    parser.add_argument('-e', '--evidence_file',
                        help='Name of tsv file with the priority score',
                        type=str, required=True)
    parser.add_argument('-c', '--cell_types_file',
                        help='Name of tsv file with cell line names per cancer type',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file. (gzip compressed json)',
                        type=str, required=True)
    parser.add_argument('-l', '--log_file',
                        help='Name of log file. If not provided logs are written to standard error.',
                        type=str, required=False)

    args = parser.parse_args()
    desc_file = args.descriptions_file
    evid_file = args.evidence_file
    cell_file = args.cell_types_file
    out_file = args.output_file

    # If no logfile is specified, logs are written to stderr 
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.log_file:
        logging.config.fileConfig(filename=args.log_file)
    else:
        logging.StreamHandler(sys.stderr)

    # Log parameters:
    logging.info(f'Evidence file: {evid_file}')
    logging.info(f'Description file: {desc_file}')
    logging.info(f'Cell type annotation: {cell_file}')
    logging.info(f'Output file: {out_file}')

    # Read files:
    evidence_df = pd.read_csv(evid_file, sep='\t')
    description_df = pd.read_csv(desc_file, sep='\t')
    cell_lines_df = pd.read_csv(cell_file, sep='\t')

    # Logging dataframe stats:
    logging.info(f'Number of evidence: {len(evidence_df)}')
    logging.info(f'Number of descriptions: {len(description_df)}')
    logging.info(f'Number of cell/tissue annotation: {len(cell_lines_df)}')

    ##
    ## Process annotation
    ##

    # Merging description with cell types and tissue:
    tissue_desc = description_df.merge(cell_lines_df[['Name','Tissue']], left_on='tissue_or_cancer_type', how='inner', right_on='Tissue')
    cell_desc = description_df.merge(cell_lines_df[['Name','Cancer Type']], left_on='tissue_or_cancer_type', how='inner', right_on='Cancer Type')

    # Concatenating annotation:
    merged_annotation = pd.concat([tissue_desc, cell_desc], ignore_index=True)

    # Aggregating names accross disease/targets:
    pooled_annotation = merged_annotation.groupby(['efo_id', 'tissue_or_cancer_type', 'method']).agg(
            {'Name': lambda x: list(x)}).reset_index()

    # Updating columns:
    pooled_annotation = (
        pooled_annotation
        .drop(['method',],axis=1)
        .rename(columns={
            'efo_id':'diseaseFromSourceMappedId',
            'Name': 'diseaseCellLines',
            'tissue_or_cancer_type': 'diseaseFromSource',
        })
    )

    ##
    ## Process evidence file:
    ##

    # Some columns from the evidence file are not needed:
    evidence_df = (
        evidence_df
        .drop(['pmid','gene_set_name', 'disease_name'], axis=1)
        .rename(columns={
            'target_id': 'targetFromSourceId', 
            'disease_id':'diseaseFromSourceMappedId',
            'score': 'resourceScore',
        })
    )
    
    # Replace some target ids: 
    evidence_df.targetFromSourceId = evidence_df.targetFromSourceId.apply(lambda x: CRISPR_SYMBOL_MAPPING[x] if x in CRISPR_SYMBOL_MAPPING else x)

    ##
    ## Annotate evidence:
    ##

    # Merging evidence and annotations:
    annotated_evidence = evidence_df.merge(pooled_annotation, on='diseaseFromSourceMappedId', how='outer',indicator=True)

    # Checking if all disease terms got matched:
    if len(annotated_evidence.loc[annotated_evidence._merge != 'both']) == 0:
        logging.info('Cell/tissue annotation and evidence successfully merged.')
    else:
        logging.warning('Problems with matching diseases between annotation and evidence. The following rows were problematic:')
        logging.warning(annotated_evidence.loc[annotated_evidence._merge != 'both'])
        annotated_evidence = annotated_evidence[annotated_evidence._merge != 'both']
        logging.warning(f'Number of evidence with mathing diseases: {len(annotated_evidence)}')

    # Remove unused column
    annotated_evidence.drop(['_merge'], inplace=True, axis=1)

    # Update efo identifier:
    annotated_evidence.diseaseFromSourceMappedId = annotated_evidence.diseaseFromSourceMappedId.str.extract('/([^/]+?)$', expand=False)

    # Adding new columns:
    annotated_evidence['datasourceId'] = 'crispr'
    annotated_evidence['datatypeId'] = 'affected_pathway'

    logging.info(f'Saving {len(annotated_evidence)} CRISPR evidence in JSON format, GZIP compressed file: {out_file}')

    ##
    ## Save evidence:
    ## 
    annotated_evidence.to_json(out_file, compression='gzip', orient='records', lines=True)


if __name__ == "__main__":
    main()

