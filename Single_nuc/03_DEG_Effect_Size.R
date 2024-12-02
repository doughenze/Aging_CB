library(Seurat)
library(ggplot2)

setwd('')
scRNA.merge.integrated <- readRDS('sn_data/Nucseq_cere_aging_full.rds')

load('mouse_structure/cerebellum/sn_data/Nucseq_cer_AgingDEGs_results_full.bin')
largest_list <- results_list_scRNA.merge.integrated

# Effect size dotplot

average_logFC <- data.frame(CellType = character(),
                            Comparison = character(),
                            AvgLogFC = numeric())

# Loop through each cell type and comparison to calculate average logFC of positively regulated genes
for (i in seq_along(largest_list)) {
  
  cell_type <- names(largest_list)[i]
  print(cell_type)
  comparisons <- names(largest_list[[i]])
  
  for (comparison in comparisons) {
    # Access the dataframe within each comparison
    ressig <- largest_list[[i]][[comparison]]$ressig
    
    # Check if the 'avg_log2FC' column exists
    if ("avg_log2FC" %in% colnames(ressig)) {
      # Filter for positively regulated DEGs
      pos_genes <- ressig[ressig$avg_log2FC > 0, ]
      
      # Order by 'avg_log2FC' and select the top 10
      top_pos_genes <- pos_genes[order(pos_genes$avg_log2FC, decreasing = TRUE), ]
      if (nrow(top_pos_genes) > 10) {
        top_pos_genes <- top_pos_genes[1:10, ]
      }
      
      # Calculate average log fold change of these top 10 genes
      avg_logFC <- mean(top_pos_genes$avg_log2FC, na.rm = TRUE)
      
      # Append to the data frame
      average_logFC <- rbind(average_logFC, data.frame(CellType = cell_type,
                                                       Comparison = comparison,
                                                       AvgLogFC = avg_logFC))
    }
  }
}

average_logFC$AvgLogFC <- ifelse(is.nan(average_logFC$AvgLogFC), 0, average_logFC$AvgLogFC)

# Convert 'Comparison' to a factor with a specific order if needed
average_logFC$Comparison <- factor(average_logFC$Comparison, levels = c("12m_vs_3m", "18m_vs_3m", "24m_vs_3m"))

ggplot(average_logFC, aes(x = Comparison, y = AvgLogFC, color = CellType)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  # Add jitter to points
  theme_minimal() +
  labs(title = "Average Log Fold Change of Top 10 Positively Regulated DEGs by Cell Type",
       x = "Age Comparison",
       y = "Average Log Fold Change",
       color = "Cell Type") +
  scale_color_manual(values = c25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))