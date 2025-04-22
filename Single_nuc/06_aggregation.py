#Author James Haberberger
import pathlib
import pandas as pd
import scanpy as sc
import anndata as ad

adata = sc.read_h5ad("/scratch/users/andytsai/M1_M16/data/raw/full_dataset.h5ad")

demographics = pd.read_excel("/scratch/users/andytsai/M1_M16/data/raw/demographics.xlsx")
adata.obs = pd.merge(adata.obs, demographics[["Microlgia ID", "Donor ID"]], left_on="specimen", right_on="Microlgia ID", how="left")

adata.write_h5ad("/scratch/users/andytsai/M1_M16/data/processed/1_aggregation.h5ad")