import numpy as np
import pandas as pd
import scanpy as sc
import os
import argparse


class PseudobulkExpression:
    """Collection of steps to process single-cell expression into pseudobulked data."""
    
    def read_h5ad(self):
        """Read and preprocess source h5ad data."""
        self.adata = sc.read(self.h5ad_source_data_path)
        # Ensembl gene IDs are used as index throughout the processing code.
        self.adata.var_names_make_unique()
        # Check AnnData object for some basic requirements.
        if self.adata.var.index.name != 'ensg':
            raise ValueError(f"Expected index name 'ensg', got '{adata.var.index.name}'")
        

    def filter_h5ad(self,min_cells,min_genes,method='10X'):
        """Filter out undesired cells and genes.
        min_cells: Minimum number of cells a gene must be expressed in.
        min_genes: Minimum number of genes a cell must express.
        method: The method used to generate the data, e.g. 10X, SmartSeq2.
        """
        adata = self.adata
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata = adata[adata.obs['method'] == method]

    def normalise_h5ad(self,method='logcp10k'):
        """Normalise the AnnData object.
        method: The normalisation method e.g. logcp10k, scTransform. If it already exists in adata.layers, it will be used.
        """
        adata = self.adata
        adata.layers['counts'] = adata.X.copy()
        if method in adata.layers:
            adata.X = adata.layers[method]
        elif method == 'logcp10k':
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    
    def combine_annotations(self,annotation_colnames):
        """Combine annotations.
        annotation_colnames: List of column names to combine.
        """
        adata = self.adata
        combined_colname = '::'.join(annotation_colnames)
        # Combine annotation columns into a single column.
        adata.obs[combined_colname] = adata.obs[annotation_colnames].astype(str).agg('::'.join, axis=1)

    def pseudobulk_data(self,aggregation_colname,donor_colname,min_cells,method='dMean'):
        """Calculate pseudobulk data.
        aggregation_colname: The annotation column name to aggregate on.
        donor_colname: The donor column name to aggregate on.
        min_cells: Minimum number of cells that an annotation-donor combination must have to be included.
        method: The method used to aggregate the data, e.g. dMean, dSum.
        """
        # Code adapted from https://github.com/wtsi-hgi/QTLight/blob/main/bin/aggregate_sc_data.py
        adata = self.adata
        print(aggregation_colname)
        print("----------")
        try:
            aggregation_column = self.adata.obs[aggregation_colname]
        except KeyError:
            print(f"Aggregation column {aggregation_colname} doesn't exist in AnnData object")
            return
        # Loop through each unique annotation in the aggregation column.
        for annotation in aggregation_column.unique():
            aggregated_data=pd.DataFrame()
        
            print("----------")
            adata.strings_to_categoricals()

            # Slice the AnnData object to only include the current annotation.
            annot_adata = adata[adata.obs[aggregation_colname]==annotation]
            annot_index = set(adata[adata.obs[aggregation_colname]==annotation].obs.index)
            aggregated_data_pre=pd.DataFrame()
            for donor in annot_adata.obs[donor_colname].unique():
                donor_index = set(adata[adata.obs[donor_colname]==donor].obs.index)
                annot_donor_index = set(annot_index.intersection(donor_index))
                donor_adata = adata[adata.obs.index.isin(annot_donor_index)]
                if donor_adata.obs.shape[0] >= min_cells:
                    if (method =='dSum'):
                        f = donor_adata.to_df()
                        data_aggregated_for_annot_and_individual = pd.DataFrame(f.sum(axis = 0))
                        data_aggregated_for_annot_and_individual.set_index(f.columns,inplace=True)
                    elif (method =='dMean'):
                        f = donor_adata.to_df()
                        data_aggregated_for_annot_and_individual = pd.DataFrame(f.mean(axis = 0))
                        data_aggregated_for_annot_and_individual.set_index(f.columns,inplace=True)
                    else:
                        raise ValueError('Wrong method specified, please use dMean or dSum')
                    data_aggregated_for_annot_and_individual.rename(columns={0:donor},inplace=True)
                    aggregated_data_pre=pd.concat([aggregated_data_pre,data_aggregated_for_annot_and_individual],axis=1)
            aggregated_data=pd.concat([aggregated_data,aggregated_data_pre],axis=1)


            output_dir = 'results/pseudobulk/'
            os.makedirs(output_dir, exist_ok=True)
            aggregated_data.to_csv(f'{output_dir}{method}-{aggregation_colname}.tsv',
                                   sep='\t', index=True)

    def main(self):
        self.read_h5ad()
        self.filter_h5ad(min_cells=3,min_genes=200)
        self.normalise_h5ad()
        self.combine_annotations(['tissue','cell_type'])
        for annotation in ['tissue::cell_type', 'tissue', 'cell_type']:
            for method in ['dMean','dSum']:
                self.pseudobulk_data(aggregation_colname=annotation,
                                     donor_colname='donor',
                                     min_cells=5,
                                     method=method)
        print("Pseudobulk expression completed")

    def __init__(self, h5ad_source_data_path):
        self.h5ad_source_data_path = h5ad_source_data_path

parser = argparse.ArgumentParser()
parser.add_argument(
    "--h5ad-source-data-path", required=True, type=str, help="Input h5ad file."
)

if __name__ == "__main__":
    args = parser.parse_args()
    PseudobulkExpression(args.h5ad_source_data_path).main()
    

