#Author James Haberberger
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

def annotate_gene_groups(adata):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    return adata



def generate_qc_plots(adata):
    sc.pl.violin(
        adata,
        [
            "total_counts", 
            "n_genes_by_counts", 
            "pct_counts_mt"
        ],
        jitter=0.4,
        multi_panel=True,
    )

    sc.pl.scatter(
        adata, 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_mt"
    )

    sc.pl.scatter(
        adata, 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_ribo"
    )
    
    sc.pl.scatter(
        adata, 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_hb"
    )

def identify_outliers(adata):

    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ribo", "hb"], 
        inplace=True,
        percent_top=[20],
        log1p=True
    )

    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    
    return adata 

adata = sc.read_h5ad("/scratch/users/andytsai/M1_M16/data/processed/1_aggregation.h5ad")

adata = annotate_gene_groups(adata)
generate_qc_plots(adata)
adata = identify_outliers(adata)

adata = adata[
    (~adata.obs.outlier) 
    & (~adata.obs.mt_outlier)
].copy()
adata.write_h5ad("/scratch/users/andytsai/M1_M16/data/processed/2_initial_preprocessing.h5ad")