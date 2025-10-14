library(Seurat)
library(ggplot2)
library(dplyr)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

setwd('/oak/stanford/groups/quake/doug/bruno_transfer/Papers/Aging_CB/Aging_CB/Single_nuc/')
scRNA.merge.integrated <- readRDS('snRNAseq_LogNorm_integrated_mergeCluster_finectyped.rds')

load('Nucseq_cer_AgingDEGs_results_full_rev.bin')

largest_list <- results_list_scRNA.merge.integrated
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
average_logFC$Comparison <- factor(average_logFC$Comparison, levels = c("12_vs_3", "18_vs_3", "24_vs_3"))

ggplot(average_logFC, aes(x = Comparison, y = AvgLogFC, color = CellType)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  # Add jitter to points
  theme_minimal() +
  labs(title = "Average Log Fold Change of Top 10 Positively Regulated DEGs by Cell Type",
       x = "Age Comparison",
       y = "Average Log Fold Change",
       color = "Cell Type") +
  scale_color_manual(values = c25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text for readability
