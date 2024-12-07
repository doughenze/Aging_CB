library(Seurat)
library(ggplot2)

setwd('')
scRNA.merge.integrated <- readRDS('snRNAseq_LogNorm_integrated_mergeCluster.rds')

scRNA.merge.integrated@meta.data$Age <- gsub("\\D", "", scRNA.merge.integrated@meta.data$orig.ident)

scRNA.merge.integrated@meta.data$cell_type_age <-
  paste(scRNA.merge.integrated@meta.data$celltype,
        scRNA.merge.integrated@meta.data$Age,
        sep = '-'
  )
DefaultAssay(scRNA.merge.integrated) <- 'RNA'
#and the default idents to our newly-set up cell_type_age 
Idents(scRNA.merge.integrated) <- 'cell_type_age'
#Test used is MAST
test_to_use <- 'MAST'
#We also define as the minimal percent a gene needs to be expressed as 5% in at least one of the compared cell groups 
min_pct <- 0.05
logfc_threshold <- 0.2
assay_to_use <- 'RNA'
latent.vars = NULL
random.seed <- 123


#We'll setup a list that will store the results from all cell_types tested
results_list_scRNA.merge.integrated <- list()
#scRNA.merge.integrated <- JoinLayers(scRNA.merge.integrated, by = "Age")

for (cell_type_to_test in unique(scRNA.merge.integrated$celltype)) {
  # We'll set up a temporary list to store the results tables
  results_list_cell_typeLevel <- list()
  print(as.character(cell_type_to_test))
  
  # We set the ages we want to compare
  age_comparisons <- list(c("24", "3"), c("18", "3"), c("12", "3"))
  
  # Loop through each age comparison
  for (ages in age_comparisons) {
    cond1 <- ages[1]
    cond2 <- ages[2]
    
    # concatenate that with the currently analysed cell_type
    cell_type_cond1 <- paste(cell_type_to_test, cond1, sep = '-')
    cell_type_cond2 <- paste(cell_type_to_test, cond2, sep = '-')
    print(paste(cell_type_cond1, 'vs', cell_type_cond2))
    
    # Perform the differential expression analysis
    diffmarkers_table <-
      FindMarkers(
        object = scRNA.merge.integrated,
        ident.1 = cell_type_cond1,
        ident.2 = cell_type_cond2,
        test.use = test_to_use,
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        assay = assay_to_use,
        verbose = TRUE,
        latent.vars = latent.vars
      )
    
    # Store the genes (currently only in rownames) to a separate column
    diffmarkers_table$gene_symbol <- row.names(diffmarkers_table)
    
    # The default p-adjust method for MAST is Bonferroni correction, we will also perform the BH correction
    diffmarkers_table$p_adjBH <- p.adjust(diffmarkers_table$p_val, method = 'BH')
    
    # Write out the comparison so that we can access the results table easily later
    comparison <- paste(cond1, cond2, sep = '_vs_')
    
    # Then we store the table containing the differential expression results in the results list
    results_list_cell_typeLevel[[comparison]]$resall <- diffmarkers_table
    
    # For easier inspection of the results, we also store a subset of the table, where we already selected for genes passing the significance threshold of 0.05
    resSig <- subset(diffmarkers_table, p_val_adj < 0.05)
    results_list_cell_typeLevel[[comparison]]$ressig <- resSig
  }
  
  # We store this temporary results list in the main list we created before starting the loop and save this before starting with the next cell_type
  results_list_scRNA.merge.integrated[[as.character(cell_type_to_test)]] <- results_list_cell_typeLevel
  save(results_list_scRNA.merge.integrated, 
       file='mouse_structure/cerebellum/sn_data/Nucseq_cer_AgingDEGs_results_full_rev.bin')
}