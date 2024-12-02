# This work flow is adapted mainly from: https://github.com/OliInTheValley/SpatioTemporal_Analysis/tree/main/Rscripts

library(Seurat)
library(ggplot2)

# set the working directory

setwd('')

# six different datsets young old male female, neuron and non-neuron
dataset_1_dir <- "run_count_CB_3FN/outs/filtered_feature_bc_matrix/"
data_young_01 <- Read10X(data.dir = dataset_1_dir)
data_young_01 <- CreateSeuratObject(counts = data_young_01, min.cells = 3, min.features = 40)
data_young_01$sampleID <- '3m_Fem'


dataset_2_dir <- "run_count_CB_3MN/outs/filtered_feature_bc_matrix/"
data_young_02 <- Read10X(data.dir = dataset_2_dir)
data_young_02 <- CreateSeuratObject(counts = data_young_02, min.cells = 3, min.features = 40)
data_young_02$sampleID <- '3m_Male'

dataset_3_dir <- "run_count_CB_3g/outs/filtered_feature_bc_matrix/"
data_young_03 <- Read10X(data.dir = dataset_3_dir)
data_young_03 <- CreateSeuratObject(counts = data_young_03, min.cells = 3, min.features = 40)
data_young_03$sampleID <- '3m_Nonneuron'

dataset_4_dir <- "run_count_CB_24FN/outs/filtered_feature_bc_matrix/"
data_old_01 <- Read10X(data.dir = dataset_4_dir)
data_old_01 <- CreateSeuratObject(counts = data_old_01, min.cells = 3, min.features = 40)
data_old_01$sampleID <- '24m_Fem'


dataset_5_dir <- "run_count_CB_24MN/outs/filtered_feature_bc_matrix/"
data_old_02 <- Read10X(data.dir = dataset_5_dir)
data_old_02 <- CreateSeuratObject(counts = data_old_02, min.cells = 3, min.features = 40)
data_old_02$sampleID <- '24m_Male'

dataset_6_dir <- "run_count_CB_24g/outs/filtered_feature_bc_matrix/"
data_old_03 <- Read10X(data.dir = dataset_6_dir)
data_old_03 <- CreateSeuratObject(counts = data_old_03, min.cells = 3, min.features = 40)
data_old_03$sampleID <- '24m_Nonneuron'


dataset_7_dir <- "run_count_CB_12FN/outs/filtered_feature_bc_matrix/"
data_mid_01 <- Read10X(data.dir = dataset_7_dir)
data_mid_01 <- CreateSeuratObject(counts = data_mid_01, min.cells = 3, min.features = 40)
data_mid_01$sampleID <- '12m_Fem'


dataset_8_dir <- "run_count_CB_12MN/outs/filtered_feature_bc_matrix/"
data_mid_02 <- Read10X(data.dir = dataset_8_dir)
data_mid_02 <- CreateSeuratObject(counts = data_mid_02, min.cells = 3, min.features = 40)
data_mid_02$sampleID <- '12m_Male'

dataset_9_dir <- "run_count_CB_12g/outs/filtered_feature_bc_matrix/"
data_mid_03 <- Read10X(data.dir = dataset_9_dir)
data_mid_03 <- CreateSeuratObject(counts = data_mid_03, min.cells = 3, min.features = 40)
data_mid_03$sampleID <- '12m_Nonneuron'

dataset_10_dir <- "run_count_CB_18FN/outs/filtered_feature_bc_matrix/"
data_18m_01 <- Read10X(data.dir = dataset_10_dir)
data_18m_01 <- CreateSeuratObject(counts = data_18m_01, min.cells = 3, min.features = 40)
data_18m_01$sampleID <- '18m_Fem'

dataset_11_dir <- "run_count_CB_18MN/outs/filtered_feature_bc_matrix/"
data_18m_02 <- Read10X(data.dir = dataset_11_dir)
data_18m_02 <- CreateSeuratObject(counts = data_18m_02, min.cells = 3, min.features = 40)
data_18m_02$sampleID <- '18m_Male'

dataset_12_dir <- "run_count_CB_18g/outs/filtered_feature_bc_matrix/"
data_18m_03 <- Read10X(data.dir = dataset_12_dir)
data_18m_03 <- CreateSeuratObject(counts = data_18m_03, min.cells = 3, min.features = 40)
data_18m_03$sampleID <- '18m_Nonneuron'


Nuc_Seq_list <- list()
Nuc_Seq_list[['3m_Fem']] <- data_young_01
Nuc_Seq_list[['3m_Male']] <- data_young_02
Nuc_Seq_list[['3m_Nonneuron']] <- data_young_03
Nuc_Seq_list[['12m_Fem']] <- data_mid_01
Nuc_Seq_list[['12m_Male']] <- data_mid_02
Nuc_Seq_list[['12m_Nonneuron']] <- data_mid_03
Nuc_Seq_list[['18m_Fem']] <- data_18m_01
Nuc_Seq_list[['18m_Male']] <- data_18m_02
Nuc_Seq_list[['18m_Nonneuron']] <- data_18m_03
Nuc_Seq_list[['24m_Fem']] <- data_old_01
Nuc_Seq_list[['24m_Male']] <- data_old_02
Nuc_Seq_list[['24m_Nonneuron']] <- data_old_03

## now we are going to loop through the list just like was done in Hahn et al

genes_to_remove <- c()

#Loop over the other datasets and place in a list
for (sample_to_analyze in names(Nuc_Seq_list) ) {
  print(sample_to_analyze)
  #We'll extract the respective single cell object
  sc_RNAseq_object <- Nuc_Seq_list[[sample_to_analyze]]
  #Here we remove any genes that we defined before starting the loop
  if (length(genes_to_remove) > 0) {
    #We extract the count matrix of the respective sample
    counts <- GetAssayData(sc_RNAseq_object, assay = "RNA")
    #subset the count matrix to all genes except those in the genes_to_remove vector
    counts <- counts[-(which(rownames(counts) %in% genes_to_remove)),]
    #Then subset the scRNA object to contain only genes that exist in the subsetted count matrix
    sc_RNAseq_object <- subset(sc_RNAseq_object,features = rownames(counts) )
  }
  #We re-create the sc object
  #cleaned_matrix <- as.matrix(sc_RNAseq_object@assays$RNA@layers$counts)
  #sc_RNAseq_object <- CreateSeuratObject(counts = cleaned_matrix, min.cells = 3, min.features = 50)
  sc_RNAseq_object@meta.data$sample <- sample_to_analyze
  #we quantify the percent of reads mapping to mitochondrial and ribosomal genes. Mitochondrial reads are an indicator of damaged cells/nuclei
  sc_RNAseq_object[["percent.mt"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "mt-")
  sc_RNAseq_object[["percent.rpl"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "rp-")
  #We filter for cells/nuclei that express more than 400 genes but less than 7000, and retain nuclei with less than 5% mitochondrial reads
  sc_RNAseq_object <- subset(sc_RNAseq_object, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 5)
  #We'll assign the sample name into a new column 'sampleID'
  sc_RNAseq_object$sampleID <- sample_to_analyze
  
  sc_RNAseq_object <- ScaleData(sc_RNAseq_object, features = rownames(sc_RNAseq_object))
  
  sc_RNAseq_object[['sampleID_sex']] <- sub(".*_(.*)", "\\1", sc_RNAseq_object$sampleID)
  
  #Now we assign the object back into the main list 
  Nuc_Seq_list[[sample_to_analyze]] <- sc_RNAseq_object
  rm(sc_RNAseq_object)
  
}


# time to integrate the data together now
scRNA.merge <- Reduce(function(x,y) merge(x,y) , Nuc_Seq_list) 
rm(Nuc_Seq_list)

# Now we'll integrate the transcriptomes to accomplish a better harmonization across all individual samples
###Integration
DefaultAssay(scRNA.merge) <- 'RNA'
#We utilize Seurat's standard workflow for integration using SCT. First we'll split the object back into its original samples and perform SCTransform on these
object_splitlist <- SplitObject(scRNA.merge, split.by = "sampleID")
for (i in names(object_splitlist)) {
  print(i)
  object_splitlist[[i]] <- SCTransform(object_splitlist[[i]], verbose = T, assay = 'RNA', vars.to.regress = c('nCount_RNA', 'nFeature_RNA'), )
}
#We utilitze the selctIntegrationFeatures (with 500 anchor genes) and PrepSCTIntegration with default settings
Integration.features <- SelectIntegrationFeatures(object.list = object_splitlist, nfeatures = 500)
object_splitlist <- PrepSCTIntegration(object.list = object_splitlist, anchor.features = Integration.features, verbose = T)
#the final Integration is then conducted
integration.anchors <- FindIntegrationAnchors(object.list = object_splitlist, normalization.method = "SCT",
                                              anchor.features = Integration.features, verbose = T)
#Having identified the anchors, we can now run the integration
scRNA.merge.integrated <- IntegrateData(anchorset =integration.anchors, normalization.method = "SCT",
)

#Before we move forward, we do a bit of relabeling and add informatino about the age groups
scRNA.merge.integrated@meta.data$Age <-  plyr::mapvalues(x = scRNA.merge.integrated@meta.data$sampleID,
                                                         from = c('3m_Fem',
                                                                  '3m_Male',
                                                                  '3m_Nonneuron',
                                                                  '12m_Fem',
                                                                  '12m_Male',
                                                                  '12m_Nonneuron',
                                                                  '18m_Fem',
                                                                  '18m_Male',
                                                                  '18m_Nonneuron',
                                                                  '24m_Fem',
                                                                  '24m_Male',
                                                                  '24m_Nonneuron'
                                                         ),
                                                         to = c('3m',
                                                                '3m',
                                                                '3m',
                                                                '12m',
                                                                '12m',
                                                                '12m',
                                                                '18m',
                                                                '18m',
                                                                '18m',
                                                                '24m',
                                                                '24m',
                                                                '24m')
)


#it is a bit easier to work with for now to use discrete values for the age groups instead of using age as a continious variable
scRNA.merge.integrated@meta.data$Age <- factor(scRNA.merge.integrated@meta.data$Age, 
                                               levels = c( 
                                                 '3m',
                                                 '12m',
                                                 '18m',
                                                 '24m'
                                               ), 
                                               ordered = T)


#we'll also extract the replicate information
scRNA.merge.integrated@meta.data$replicate <- unlist(lapply(strsplit((scRNA.merge.integrated@meta.data$sampleID), split = "_"), function(x) x[[2]]))
#and the cell BC
scRNA.merge.integrated@meta.data$index <- unlist(lapply(strsplit(row.names(scRNA.merge.integrated@meta.data), split = "-"), function(x) x[[1]]))
#Finlally, we'll create a composite cell index based on sampleID and cell BC
scRNA.merge.integrated@meta.data$cell_index <- paste(scRNA.merge.integrated@meta.data$tissue, scRNA.merge.integrated@meta.data$Age, scRNA.merge.integrated@meta.data$replicate, scRNA.merge.integrated@meta.data$index, sep = '_')

#We'll also reassing some information that we already have to make the data better match to the metadata of bulk and spatial data
scRNA.merge.integrated@meta.data$tissue <- "Cerebellum"
scRNA.merge.integrated@meta.data$sex <- scRNA.merge.integrated@meta.data$sampleID_sex


#Now we will run the dimensionality reduction and clustering on the integrated data slot
DefaultAssay(scRNA.merge.integrated)
scRNA.merge.integrated <- RunPCA(object = scRNA.merge.integrated, verbose = T)
scRNA.merge.integrated <- FindNeighbors(scRNA.merge.integrated, dims = 1:15)
scRNA.merge.integrated <- FindClusters(scRNA.merge.integrated, resolution = 1)
scRNA.merge.integrated <- RunUMAP(scRNA.merge.integrated, dims = 1:15, verbose = T)

#We can also look just at transcriptome as a umap and visualize if and how similar anterior/posterior transcriptomes are
DimPlot(scRNA.merge.integrated, label = T)
DimPlot(scRNA.merge.integrated, group.by = 'sampleID')

#Now we'll run SCTransformation as a way to get the expression data normalized in a way that we can analyze gene expression, calculate signature scores etc.
saveRDS(scRNA.merge.integrated, file = 'sn_data/Nucseq_cere_aging_full.rds')
scRNA.merge.integrated <- readRDS('sn_data/Nucseq_cere_aging_full.rds')

#First we need to re-set the default assay ("integrated" is only sufficient for the dimensionality reduction)
DefaultAssay(scRNA.merge.integrated) <- "RNA"
#Run SCT on the whole dataset
scRNA.merge.integrated <- SCTransform(scRNA.merge.integrated, assay = "RNA", verbose = T)
#In addtion, we'll also run the standard data normalization and scaling
DefaultAssay(scRNA.merge.integrated) <- "RNA"

scRNA.merge.integrated <- NormalizeData(scRNA.merge.integrated)
scRNA.merge.integrated <- FindVariableFeatures(scRNA.merge.integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA.merge.integrated)
scRNA.merge.integrated <- ScaleData(scRNA.merge.integrated, features = all.genes)

#Let's save the integrated dataset
saveRDS(scRNA.merge.integrated, file = 'sn_data/Nucseq_cere_aging_full.rds')

scRNA.merge.integrated <- readRDS('sn_data/Nucseq_cere_aging_full.rds')


## annotation in R
Idents(scRNA.merge.integrated) <- 'seurat_clusters'
scRNA.merge.integrated <- PrepSCTFindMarkers(object = scRNA.merge.integrated)
#MarkerList_to_annotate <- FindAllMarkers(object = scRNA.merge.integrated, only.pos = TRUE, min.pct = 0.15, 
#                                        thresh.use = 0.15, verbose = T, assay = 'SCT')

genes_of_interest <- c("Gabra6", "Ppp1r17","Eomes", "Klhl1", "Lypd6","Prkcd","Lgi2","Gdf10","Aqp4","Mobp"
                       , "Ppfibp1", "Dcn", "Kcnj8", "Ttr", "Mrc1", "C1qa", "Flt1","Foxj1") # Replace with your genes

# Assuming 'seurat_object' is your Seurat object
# Create the dot plot
dotplot <- DotPlot(scRNA.merge.integrated, features = genes_of_interest) + 
  scale_color_gradientn(colors = c("blue", "white", "red")) + # Optional: Adjust color gradient
  theme_minimal() + 
  ggtitle("Expression of Genes Across Clusters")

# Print the dot plot
print(dotplot)

DimPlot(scRNA.merge.integrated, label = T)

current.cluster.ids <- c(0, 1, 2, 3,
                         4, 5, 6, 7,
                         8, 9, 10, 11,
                         12, 13, 14, 15,
                         16, 17, 18, 19,
                         20, 21, 22, 23,
                         24, 25, 26, 27,
                         28)

new.cluster.ids <- c('Granule', 'Granule', 'Granule', 'ODC',
                     'ODC', 'MLI', 'PLI', 'ODC',
                     'Unknown', 'Bergmann', 'Astrocyte', 'Unknown',
                     'Microglia', 'Granule', 'Unknown', 'Unknown',
                     'Fibroblast', 'OPC', 'Unknown', 'Unknown',
                     'Unknown', 'UBC', 'Unknown', 'Purkinje',
                     'Unknown','Endothelial','Unknown', 'Choroid',
                     'Unknown')



#we use plyr's mapvalues function to re-label the clusters, 
scRNA.merge.integrated@meta.data$celltype <- plyr::mapvalues(x = scRNA.merge.integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#Let's set the default 'idents' to our cluster-level annotation
Idents(scRNA.merge.integrated) <- 'celltype'

genes_of_interest <- c("Gabra6", "Ppp1r17","Eomes", "Klhl1", "Lypd6","Prkcd","Lgi2","Gdf10","Aqp4","Mobp"
                       , "Ppfibp1", "Dcn", "Kcnj8", "Ttr", "Mrc1", "C1qa", "Flt1","Foxj1") # Replace with your genes


genes_of_interest <- c("Gabra6", "Mobp","Lypd6", "Klhl1", "Aqp4","Gdf10","Ttr","Ppfibp1","Dcn","C1qa"
                       , "Ppp1r17", "Flt1", "Kcnj8", "Eomes", "Lgi2") # Replace with your genes
genes_of_interest <- c("Ttr", "Kcnj8","Ppp1r17", "Eomes", "Ppfibp1","Dcn","C1qa","Aqp4","Gdf10","Klhl1"
                       , "Lypd6", "Mobp", "Gabra6")


dotplot <- DotPlot(scRNA.merge.integrated, features = genes_of_interest) + 
  scale_color_gradientn(colors = c("blue", "white", "red")) + # Optional: Adjust color gradient
  theme_minimal() + 
  ggtitle("Expression of Genes Across Clusters")

# Print the dot plot
print(dotplot)

#We'll remove the doublets and unknown clusters
scRNA.merge.integrated <- subset(scRNA.merge.integrated, celltype != 'Unknown')
#We'll save the object
saveRDS(scRNA.merge.integrated, file = 'sn_data/Nucseq_cere_aging_full.rds')