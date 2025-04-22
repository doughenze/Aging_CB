#Author James Haberberger
import pathlib
import pandas as pd
import scanpy as sc
import anndata as ad

adata = sc.read_h5ad("/scratch/users/andytsai/M1_M16/data/processed/2_initial_preprocessing.h5ad")

# Required to make index unique before running scrublet.
adata.obs_names_make_unique()

sc.pp.scrublet(adata, batch_key="specimen")

adata.write_h5ad("/scratch/users/andytsai/M1_M16/data/processed/3_scrublet.h5ad")
adata[adata.obs["predicted_doublet"].eq(False)].write_h5ad("/scratch/users/andytsai/M1_M16/data/processed/3_scrublet_subset.h5ad")