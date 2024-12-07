# This work flow is adapted mainly from: https://github.com/OliInTheValley/SpatioTemporal_Analysis/tree/main/Rscripts

library(Seurat)
library(ggplot2)

setwd('')

scRNA.merge.integrated <- readRDS('snRNAseq_LogNorm_integrated_mergeCluster.rds')
DimPlot(scRNA.merge.integrated)
DefaultAssay(scRNA.merge.integrated) <- "RNA"
genes_of_interest <- c("Gabra6", "Ppp1r17","Eomes", "Klhl1", "Lypd6","Prkcd","Lgi2","Gdf10","Aqp4","Mobp"
                       , "Ppfibp1", "Dcn", "Kcnj8", "Ttr", "Mrc1", "P2ry12", "Flt1","Foxj1") # Replace with your genes

# Assuming 'seurat_object' is your Seurat object
# Create the dot plot
dotplot <- DotPlot(scRNA.merge.integrated, features = genes_of_interest) + 
  scale_color_gradientn(colors = c("blue", "white", "red")) + # Optional: Adjust color gradient
  theme_minimal() + 
  ggtitle("Expression of Genes Across Clusters")

# Print the dot plot
print(dotplot)


current.cluster.ids <- c(0, 1, 3,
                         4, 5, 6, 7,
                         8, 9, 10, 11,
                         12, 13, 14, 15)

new.cluster.ids <- c('Granule', 'Oligo', 'Astrocyte',
                     'Fibroblast', 'MLI', 'Bergmann', 'Choroid',
                     'Endo', 'Endo. Mural', 'Astrocyte', 'Microglia',
                     'OPC', 'PLI','Purkinje', 'UBC')

scRNA.merge.integrated@meta.data$celltype <- plyr::mapvalues(x = scRNA.merge.integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#Let's set the default 'idents' to our cluster-level annotation
Idents(scRNA.merge.integrated) <- 'celltype'
DimPlot(scRNA.merge.integrated)
saveRDS(scRNA.merge.integrated, file = 'snRNAseq_LogNorm_integrated_mergeCluster.rds')
