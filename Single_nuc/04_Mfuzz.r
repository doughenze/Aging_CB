setwd('/home/cd973/Downloads')
library("Seurat")
library("SeuratWrappers")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("grid")
library("monocle3")
library("e1071")
library("Mfuzz")

myObj= "snRNAseq_LogNorm_integrated_mergeCluster.rds"

seuratObj<- readRDS(myObj)
seuratObj$msex<- factor(ifelse(grepl('MN',seuratObj@meta.data$orig.ident), 1, ifelse(grepl('FN',seuratObj@meta.data$orig.ident), 0, NA)))
seuratObj$age<- factor(gsub("MN|FN|G","", data.frame(do.call('rbind', strsplit(as.character(seuratObj@meta.data$orig.ident),'_',fixed=TRUE)))[,2]))
DefaultAssay(seuratObj) <- "RNA"
# create new seurat object
# seuratObj2<- CreateSeuratObject(seuratObj@assays$RNA@counts, project = "SeuratProject", assay = "RNA",  meta.data = seuratObj@meta.data)

clusters<- sort(unique(seuratObj@meta.data$seurat_clusters))
clusters<- clusters[as.numeric(as.character(clusters)) >8 ]
for(cluster in clusters){
  #cluster = 1
  print(paste0("Working with cluster-", cluster))
  subObj_1<- subset(seuratObj, subset = seurat_clusters==cluster & msex %in% c(0,1))
  pseudoBulk<- AverageExpression(object=subObj_1, group.by = "orig.ident", slot = "data")$RNA # get avg exprs
  
  expression_data <- GetAssayData(object = subObj_1, slot = "counts")
  celldata<- as.data.frame(subObj_1@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","msex","age")])
  genedata<- as.data.frame(x = row.names(subObj_1), row.names = row.names(subObj_1)); colnames(genedata) <- "gene_short_name";
  # Make the CDS object
  cds_1<- new_cell_data_set(expression_data, cell_metadata=celldata, gene_metadata=genedata)
  # cds_1<- as.cell_data_set(subObj_1)
  
  # https://cole-trapnell-lab.github.io/monocle3/docs/differential/
  # The quasi-poisson distribution doesn't have a real likelihood function.
  # The columns in results tables from evaluate_fits() and compare_models() will be NA.
  full_model<- fit_models(cds_1, model_formula_str = "~age + msex", expression_family="negbinomial")
  reduced_model<- fit_models(cds_1, model_formula_str = "~msex", expression_family="negbinomial")
  lrt_result<- compare_models(full_model, reduced_model); lrt_result$qvalue<- NULL;
  lrt_result$fdr<- p.adjust(lrt_result$p_value, method = "BH")
  write.table(lrt_result, paste('cluster', cluster, "_LRT_result.tsv", sep=""), row.names=F, quote=F, sep="\t")
  
  sig_gene1<- lrt_result[which(lrt_result$fdr<0.05),'gene_short_name']
  pseudoBulk2<- pseudoBulk[sig_gene1$gene_short_name,] # only pick significant genes from DESeq2 analysis
  geneData<- merge(unique(celldata[,c("orig.ident","age")]),t(pseudoBulk2),by.x='orig.ident',by.y='row.names',all=F)
  geneData$orig.ident<- NULL;
  geneData<- aggregate(.~ age, data=geneData, FUN = mean)
  geneData<- as.data.frame(t(data.frame(geneData,row.names=1)))
  geneData<- geneData[,c("3","12","18","24")]
  
  exprs<- function(eset){Biobase::exprs(eset)} #mfuzz errors
  eset <- new("ExpressionSet",exprs = as.matrix(geneData)) #Biobase function
  eset <- standardise(eset)
  m <- mestimate(eset) # get optimal m value for Mfuzz

  if(length(sig_gene1$gene_short_name)>=9){
    clust <- mfuzz(eset, c = 9, m = m) #clustering
    }else{ # if gene number < 9, can nnot split into 9 clusters
      clust <- mfuzz(eset, c = length(sig_gene1$gene_short_name), m = m) #clustering
    }

  library(RColorBrewer)
  color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
  pdf(paste('cluster', cluster, '_Mfuzz_vis.pdf', sep=''))
  mfuzz.plot(eset, clust,mfrow=c(3,3),new.window= FALSE,time.labels=colnames(geneData),colo = color.2)
  dev.off()
  
  clust_data<- data.frame(Geneid= names(clust$cluster), mfuzz_cluster=as.numeric(clust$cluster))
  clust_data2<- merge(lrt_result, clust_data,by.x='gene_short_name',by.y='Geneid', all.x=T)
  write.csv(clust_data2, paste('cluster', cluster, '_Mfuzz_cluster_result.csv', sep=''), row.names=F, quote=F)
  
}