library(Seurat)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

setwd('/oak/stanford/groups/quake/doug/bruno_transfer/Papers/Aging_CB/Aging_CB/Single_nuc/')
# ------------------------------------------------------------------
# Load the integrated Seurat object
# ------------------------------------------------------------------
seurat_obj <- readRDS("snRNAseq_LogNorm_integrated_mergeCluster.rds")

# Remove the doublets

DefaultAssay(seurat_obj) <- "RNA"

# ------------------------------------------------------------------
# 1.  DOUBLETS WITH Doublet Finder
# ------------------------------------------------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# 2. Parameter sweep to pick an optimal pK
sweep.res  <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats<- summarizeSweep(sweep.res, GT = FALSE)
pK         <- find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]

# 3. Estimate expected doublet rate
nExp <- round(0.08 * ncol(seurat_obj))        # 5 % as a starting point
homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)
nExp.adj <- round(nExp * (1 - homotypic.prop))

# 4. Run DoubletFinder
seurat_obj <- doubletFinder(
  seurat_obj,
  PCs            = 1:20,
  pN             = 0.25,
  pK             = 0.09,
  nExp           = nExp.adj,#.adj,
  #reuse.pANN     = FALSE,
  sct            = FALSE)

# 5. Add the call to metadata
df_col <- grep("DF.classification", colnames(seurat_obj@meta.data), value = TRUE)
seurat_obj$DF_doublet <- seurat_obj@meta.data[[df_col]]

DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "DF_doublet",
  pt.size = 0.1
)

seurat_obj <- subset(seurat_obj, subset = DF_doublet == "Singlet")

microglia <- subset(seurat_obj, subset = celltype == "Microglia")

microglia <- FindVariableFeatures(microglia)                      # NEW
microglia <- ScaleData(microglia)                                 # NEW
microglia <- RunPCA(microglia)                                    # NEW
microglia <- FindNeighbors(microglia, dims = 1:20)                # NEW
microglia <- FindClusters(microglia, resolution = 0.3)            # NEW
microglia <- RunUMAP(microglia, dims = 1:20)                      # NEW

current.cluster.ids <- c(0, 1, 2, 3, 4)

microglia <- c('Microglia', 'Microglia', 'BAM','Oligo', 'Microglia')

microglia@meta.data$subtype <- plyr::mapvalues(x = microglia@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#Let's set the default 'idents' to our cluster-level annotation
Idents(microglia) <- 'subtype'

seurat_obj$micro_subtype <- NA
seurat_obj$micro_subtype[Cells(microglia)] <- as.character(microglia$subtype)

# --------------------------------------------------------------
# 3.  DEFINE A REFINED CELL-TYPE COLUMN            ---- CHANGED
# --------------------------------------------------------------
# use microglia subtypes where present, else keep original `celltype`
seurat_obj$celltype_refined <- ifelse(
  is.na(seurat_obj$micro_subtype),
  as.character(seurat_obj$celltype),
  as.character(seurat_obj$micro_subtype)
)

# --------------------------------------------------------------
# 4.  SET ACTIVE IDENTITIES TO THE REFINED COLUMN   ---- CHANGED
# --------------------------------------------------------------
Idents(seurat_obj) <- "celltype_refined"

saveRDS(seurat_obj,'snRNAseq_LogNorm_integrated_mergeCluster_finectyped.rds')










####
#GENERATE THE QC IMAGES AND THE NEW UMAP
####

seurat_obj <- readRDS("snRNAseq_LogNorm_integrated_mergeCluster_finectyped.rds")

####
# Figure 1a
####
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

user_pal <- c(
  "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78",
  "#2CA02C", "#98DF8A", "#D62728", "#FF9896",
  "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
  "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7",
  "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5",
  "#7F9EC0", "#FFAB60", "#5AB4AC"
)

# -----------------------------------------------------------
# 1.  PREP DATA
# -----------------------------------------------------------
ctype   <- Idents(seurat_obj) |> as.character()
umap_df <- Embeddings(seurat_obj, "umap") |>
  as.data.frame() |>
  setNames(c("UMAP1", "UMAP2")) |>
  transform(celltype = ctype)

ctype_counts <- umap_df |>
  dplyr::count(celltype, name = "nuclei") |>
  dplyr::arrange(desc(nuclei)) |>
  transform(label = sprintf("%s  (%.2f)", celltype, nuclei / 1000))

# -----------------------------------------------------------
# 2.  ASSIGN COLOURS
# -----------------------------------------------------------
n_ct <- nrow(ctype_counts)
if (n_ct > length(user_pal)) {
  extra <- grDevices::colorRampPalette(user_pal)(n_ct)   # extend palette
  pal   <- extra[seq_len(n_ct)]
} else {
  pal <- user_pal[seq_len(n_ct)]
}
names(pal) <- ctype_counts$celltype

# -----------------------------------------------------------
# 3.  MAIN UMAP PANEL
# -----------------------------------------------------------
centroids <- umap_df |>
  dplyr::group_by(celltype) |>
  dplyr::summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, colour = celltype)) +
  ggplot2::geom_point(size = 0.4) +
  ggplot2::scale_colour_manual(values = pal, guide = "none") +
  ggrepel::geom_text_repel(
    data        = centroids,
    ggplot2::aes(label = celltype),
    colour      = "black",
    size        = 5,
    fontface    = "bold",
    segment.size= 0.2,
    max.overlaps= Inf
  ) +
  ggplot2::annotate(
    "segment", x = -20, xend = -5, y = -18, yend = -18,
    arrow = ggplot2::arrow(type = "closed", length = grid::unit(4, "pt")), size = 1
  ) +
  ggplot2::annotate(
    "segment", x = -20, xend = -20, y = -18, yend = -5,
    arrow = ggplot2::arrow(type = "closed", length = grid::unit(4, "pt")), size = 1
  ) +
  ggplot2::annotate("text", x = -4,  y = -18, label = "UMAP1",
                    hjust = 0, vjust = 1.2, size = 5, fontface = "bold") +
  ggplot2::annotate("text", x = -20.4, y = -4,  label = "UMAP2",
                    angle = 90, hjust = 0, vjust = 1.2, size = 5, fontface = "bold") +
  ggplot2::theme_void() +
  ggplot2::coord_fixed() +
  ggplot2::annotate(
    "text",
    x = min(umap_df$UMAP1) + 5,
    y = min(umap_df$UMAP2) + 5,
    label = sprintf("Nuclei: %s", scales::comma(nrow(umap_df))),
    hjust = 0, vjust = 0, size = 6, fontface = "bold"
  )

# -----------------------------------------------------------
# 4.  LEGEND PANEL
# -----------------------------------------------------------
legend_df <- transform(ctype_counts, y = rev(seq_len(n_ct)))

p_legend <- ggplot2::ggplot(legend_df, ggplot2::aes(x = 0, y = y)) +
  ggplot2::geom_point(ggplot2::aes(colour = celltype), size = 6) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    hjust = 0, nudge_x = 0.4, size = 4
  ) +
  ggplot2::scale_colour_manual(values = pal, guide = "none") +
  ggplot2::xlim(0, 5) +
  ggplot2::theme_void() +
  ggplot2::annotate(
    "text", x = 0, y = max(legend_df$y) + 1,
    label = "Nuclei count (x1000)",
    hjust = 0, vjust = 0, size = 5, fontface = "bold"
  )

# -----------------------------------------------------------
# 5. COMBINE PANELS
# -----------------------------------------------------------
patchwork::wrap_plots(p_umap, p_legend, widths = c(3, 1))

####
# Figure 1b
####
genes_of_interest <- c(
  "Aqp4", "Mrc1","Gdf10","Ttr","Pecam1","Kcnj8","Col1a1","Lgi2","Gabra6","P2ry12",
  "Lypd6","Mobp","Pdgfra","Gad2","Ppp1r17","Eomes"
)

library(scales)          # for rescale / gradient_n_pal if you want

custom_pal <- colorRampPalette(c(
  "#08306B",   # deep navy-blue
  "#2171B5",   # mid blue
  "#6BAED6",   # light blue
  "#C6DBEF",   # very light blue
  "white",     # centre
  "#FEE0D2",   # very light orange
  "#FC8D59",   # light orange-red
  "#E34A33",   # mid red
  "#B30000"    # deep crimson
))

# ------------------------------------------------------------------
# Apply to the dot-plot
# ------------------------------------------------------------------
dotplot <- DotPlot(
  seurat_obj,
  features  = genes_of_interest,
  group.by  = "celltype_refined",
  dot.scale = 8
) +
  scale_color_gradientn(
    colours = custom_pal(100),   # 100-step smooth gradient
    #limits  = c(min(dotplot$data$avg.exp.scaled),
    #            max(dotplot$data$avg.exp.scaled)),
    oob     = scales::squish      # keep values inside colour bar
  ) +
  labs(x = NULL, y = NULL, colour = "expression", size = "proportion") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(dotplot)

####


# ------------------------------------------------------------------
# Specify metadata column names (edit here if yours differ)
# ------------------------------------------------------------------
seurat_obj@meta.data$age <- gsub("\\D", "", seurat_obj@meta.data$orig.ident)
batch_col <- "orig.ident"        # e.g. "orig.ident" or another column name
age_col   <- "age"          # age metadata column (numeric or factor)

stopifnot(batch_col %in% colnames(seurat_obj@meta.data))
stopifnot(age_col   %in% colnames(seurat_obj@meta.data))

# ------------------------------------------------------------------
# QC violin plots: nFeature_RNA (genes) and nCount_RNA (UMIs) per batch
# ------------------------------------------------------------------
vln_genes <- VlnPlot(
  seurat_obj,
  features = "nFeature_RNA",
  group.by = batch_col,
  pt.size = 0.0
) + ggtitle("Genes per cell by batch") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

vln_umis <- VlnPlot(
  seurat_obj,
  features = "nCount_RNA",
  group.by = batch_col,
  pt.size = 0.0
) + ggtitle("UMIs per cell by batch") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display side-by-side
(vln_genes | vln_umis) + plot_layout(guides = "collect")

# ------------------------------------------------------------------
# UMAPs colored by batch and by age
# ------------------------------------------------------------------
umap_batch <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = batch_col,
  pt.size = 0.1
) + ggtitle("UMAP: colored by batch")

umap_age <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = age_col,
  pt.size = 0.1
) + ggtitle("UMAP: colored by age")

# Display side-by-side
umap_batch + umap_age + plot_layout(ncol = 2)


####
# Figure 1E&F
####
#install.packages("pheatmap")
library(pheatmap)

genes_to_show <- c(
  "Lilrb4a", "Lilr4b", "Zc3hav1", "Stat1", "Atp6v0d2", "Ifi207", "Cdkn2a", "Ddx60",
  "Trim30a", "Tnfsf8", "Snx20", "Cd180", "Gm4951"
)

# custom palette (same B-W-R as dot-plot)
custom_pal <- colorRampPalette(c(
  "#08306B","#2171B5","#6BAED6","#C6DBEF",
  "white",
  "#FEE0D2","#FC8D59","#E34A33","#B30000"
))

# helper to Z-score rows
scale_rows <- function(m) {
  t(scale(t(m)))
}

# ------------------------------------------------------------------
# 1.  AGE-LEVEL HEAT-MAP (all cells pooled by age)               ---
# ------------------------------------------------------------------

microglia <- subset(seurat_obj, subset = celltype_refined == "Microglia")
age_mat <- AverageExpression(
  microglia,
  features = genes_to_show,
  group.by = "age",              # <- metadata column
  layer     = "data",             # log-normalised expression
  assays   = "RNA"
)$RNA

age_mat <- scale_rows(age_mat)   # centre & scale each gene
age_mat <- age_mat[genes_to_show, ]  # enforce ordering

#row_anno <- data.frame(Cluster = gene_clusters)[genes_to_show,,drop=FALSE]

pheatmap(
  age_mat,
  color            = custom_pal(100),
  cluster_rows     = TRUE,
  cluster_cols     = FALSE,         # keep chronological order of ages
  #annotation_row   = row_anno,
  fontsize_row     = 9,
  fontsize_col     = 10,
  breaks           = seq(-3, 3, length.out = 101),
  main             = "Gene expression across ages"
)

# ------------------------------------------------------------------
# 2.  CELLTYPE-LEVEL HEAT-MAP (refined identities)               ---
# ------------------------------------------------------------------
ct_mat <- AverageExpression(
  seurat_obj,
  features = genes_to_show,
  group.by = "celltype_refined",
  layer    = "data",
  assays   = "RNA"
)$RNA

ct_mat <- scale_rows(ct_mat)
ct_mat <- ct_mat[genes_to_show, ]                     # same gene order

# optionally reorder cell types (here left-to-right as in UMAP legend)
#ct_mat <- ct_mat[, levels(seurat_obj$celltype_refined)]

pheatmap(
  ct_mat,
  color            = custom_pal(100),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,          # dendrogram for cell types
  #annotation_row   = row_anno,
  fontsize_row     = 9,
  fontsize_col     = 9,
  breaks           = seq(-3, 3, length.out = 101),
  main             = "Gene expression across refined cell types"
)




# ------------------------------------------------------------------
# 1. SUB-CLUSTER MICROGLIA
# ------------------------------------------------------------------
# --- NEW: isolate Microglia and re-run the standard Seurat workflow
microglia <- subset(seurat_obj, subset = celltype_refined == "Microglia")

microglia <- FindVariableFeatures(microglia)                      # NEW
microglia <- ScaleData(microglia)                                 # NEW
microglia <- RunPCA(microglia)                                    # NEW
microglia <- FindNeighbors(microglia, dims = 1:20)                # NEW
microglia <- FindClusters(microglia, resolution = 0.3)            # NEW
microglia <- RunUMAP(microglia, dims = 1:20)                      # NEW

# ------------------------------------------------------------------
# 2. ADD MODULE SCORES (DAM, LDAM, Homeostatic)
# ------------------------------------------------------------------
# --- NEW: define gene sets (edit as needed)
DAM_genes   <- c("Trem2", "Apoe", "Tyrobp", "Itgax", "Clec7a", "Lpl", "Cst7",
                 "Spp1", "Axl", "Cd9")
LDAM_genes  <- c("Slc25a5", "Npl", "Angptl7", "Pde2a", "Ldhb", "Cd63", "Sepp1",
                 "Sdcbp", "Adipor1", "Rbbp4", "Cndp2", "Hsd17b4", "Gpd11",
                 "Dazap2", "Hnmpk", "Rapsn", "Cat", "Kl", "Nampt", "Acsl1",
                 "Dpyd", "Cd163")
ADEM_genes <- c("Hfe", "Cyr61", "Ccl24", "Tlr7", "Cxcr6","Tap1",
                "Gcnt2","Tnf","Axl","St3gal4","Jchain","Apod",
                "Cd74", "Lpl","H2-K1", "Ltb", "Arid5a","Irak2")
Homeo_genes <- c("P2ry12","Tmem119","Cx3cr1","Sall1","Siglech","Gpr34")

microglia <- AddModuleScore(microglia, features = list(DAM_genes),   name = "DAM")
microglia <- AddModuleScore(microglia, features = list(LDAM_genes),  name = "LDAM")
microglia <- AddModuleScore(microglia, features = list(Homeo_genes), name = "Homeo")
microglia <- AddModuleScore(microglia, features = list(ADEM_genes), name = "ADEM")

# ------------------------------------------------------------------
# 3. VISUALIZE UMAP WITH MODULE SCORES
# ------------------------------------------------------------------
# --- NEW: FeaturePlot for module scores
p_dam_alt   <- FeaturePlot(microglia, "DAM1", cols = c("lightgrey","red")) + 
  ggtitle("DAM score") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"))

p_ldam_alt  <- FeaturePlot(microglia, "LDAM1", cols = c("lightgrey","red")) + 
  ggtitle("LDAM score") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"))

p_adem_alt  <- FeaturePlot(microglia, "ADEM1", cols = c("lightgrey","red")) + 
  ggtitle("ADEM score") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"))

p_homeo_alt <- FeaturePlot(microglia, "Homeo1", cols = c("lightgrey","red")) + 
  ggtitle("Homeostatic score") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"))

# Arrange without collecting guides (each plot keeps its own legend)
(p_dam_alt | p_ldam_alt | p_adem_alt | p_homeo_alt)

library(ggplot2)
(p_dam | p_ldam | p_adem | p_homeo) + plot_layout(guides = "collect")                 # NEW

# ------------------------------------------------------------------
# 4. VIOLIN PLOTS BY AGE FOR SELECTED GENES
# ------------------------------------------------------------------
# --- NEW: violin plots split by age
genes_to_plot <- c("Arhgap5", "Grin2c", "Adcy1")                              # NEW

vln_list <- lapply(
  genes_to_plot,
  function(g) VlnPlot(
    microglia,
    features = g,
    group.by = age_col,                 # defined earlier
    pt.size  = 0.1
  ) + ggtitle(g) + theme_minimal()
)

wrap_plots(vln_list)         

meta_df <- microglia@meta.data
meta_df$barcode <- rownames(meta_df)            # keep barcodes explicit

# write to disk
write.csv(
  meta_df,
  file = "microglia_metadata.csv",
  row.names = FALSE,          # barcodes now in their own column
  quote = TRUE
)
