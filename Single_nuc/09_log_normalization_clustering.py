#Author James Haberberger
import scanpy as sc

adata = sc.read_h5ad("/scratch/users/andytsai/M1_M16/data/processed/3_scrublet_subset.h5ad")

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)

# Logarithmize the data
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="specimen")
sc.pl.highly_variable_genes(adata)

# PCA
sc.tl.pca(adata)

# Harmony Integrate
sc.external.pp.harmony_integrate(adata, key="specimen")

# Neighbor Calculation
sc.pp.neighbors(adata, use_rep="X_pca_harmony")

# UMAP Calculation
sc.tl.umap(adata)

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
for res in [x*0.05 for x in range(1, 21)]:
    sc.tl.leiden(
        adata, 
        key_added=f"leiden_res_{res:4.2f}", 
        resolution=res, 
        flavor="igraph",
        n_iterations=2
    )

    sc.pl.umap(
        adata,
        color=[f"leiden_res_{res:.2f}"],
        legend_loc="on data",
    )

sc.pl.umap(
    adata,
    color=["specimen"],
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

adata.write_h5ad("/scratch/users/andytsai/M1_M16/data/processed/4_lognormalization_and_clustering.h5ad")