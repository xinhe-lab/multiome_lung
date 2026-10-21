#!/software/python-anaconda-2020.02-el7-x86_64/bin/python
"""
A python script to call celltypist to make cell type predictions
"""
#load raw count matrix and reference model
#import polars as pl
import celltypist
from celltypist import models
import scanpy as sc
import argparse
parser=argparse.ArgumentParser(description='run celltypist tools')
parser.add_argument('h5ad', type=str, help='a h5ad file for RNA raw counts matrix')
parser.add_argument('cluster_col', type=str, help='the name of the column for clustering')
parser.add_argument('out', type=str, help='path for the output table')

args=parser.parse_args()

#Show all available models that can be downloaded and used.
models.models_description()
#Download a specific model, for example, `Immune_All_Low.pkl`.
models.download_models(model = 'Immune_All_Low.pkl')


input_file=sc.read(args.h5ad) #36601x46906
#In case your input file is a gene-by-cell table.
predictions = celltypist.annotate(
  input_file, model = 'Immune_All_Low.pkl', 
  majority_voting = True,
  over_clustering = args.cluster_col,
  transpose_input = True)
#Summary information for the prediction result.
output = predictions.predicted_labels
output['conf_score'] = predictions.to_adata().obs['conf_score']
output.to_csv(args.out)
# #Examine the predicted cell type labels.
# predictions.predicted_labels
# #Examine the matrix representing the decision score of each cell belonging to a given cell type.
# predictions.decision_matrix
# #Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
# predictions.probability_matrix
# Export the three results to csv tables.
# predictions.to_table(folder = '/path/to/a/folder', prefix = '')
# Alternatively, export the three results to a single Excel table (.xlsx).
# predictions.to_table(folder = '/path/to/a/folder', prefix = '', xlsx = True)
