library(Seurat)
library(ggplot2)
library(dplyr)

setwd('/oak/stanford/groups/quake/doug/bruno_transfer/Papers/Aging_CB/Aging_CB/Single_nuc/')
seurat_obj <- readRDS('snRNAseq_LogNorm_integrated_mergeCluster_finectyped.rds')

# 1. Filter for microglia cells
mic <- subset(seurat_obj, celltype_refined == "Microglia")

# 2. Define score genes
score_genes <- c('Cxcl2','Slamf9','Cdk2', 'Atp2a3', 'Arhgap5')

# 3. Add module score (equivalent to sc.tl.score_genes)
mic <- AddModuleScore(mic, 
                      features = list(score_genes), 
                      name = "prox_score")

# Note: AddModuleScore adds a "1" suffix, so the score will be "prox_score1"
score_col <- "prox_score1"

lower_threshold <- 0
# 5. Subset to values above 0
CB_gran_assoc <- subset(mic, subset = prox_score1 >= lower_threshold)

# Optional: also get bottom 50% (equivalent to CB_far)
CB_far <- subset(mic, subset = prox_score1 < lower_threshold)

# 6. Plot score distribution (equivalent to the plotting in your function)
pdf("full_dataset_mouse_scores.pdf")
hist(scores, bins = 50, col = "lightblue", 
     main = "Prox Score Distribution", 
     xlab = "Prox Score", ylab = "Frequency")
abline(v = lower_threshold, col = "red", lty = 2, lwd = 2)
legend("topright", c("threshold"), 
       col = c("red"), lty = 2)
dev.off()

CB_gran_assoc@meta.data$Age <- gsub("\\D", "", CB_gran_assoc@meta.data$orig.ident)

# 7. Perform differential expression between young (3m) and old (24m) in top 20%
# First, subset to only 3m and 24m samples
CB_gran_assoc_filtered <- subset(CB_gran_assoc, Age %in% c("3", "24"))

# Set the identity to Age for differential expression
Idents(CB_gran_assoc_filtered) <- CB_gran_assoc_filtered$Age

# Perform differential expression (24m vs 3m)
de_results <- FindMarkers(CB_gran_assoc_filtered, 
                          ident.1 = "24", 
                          ident.2 = "3",
                          test.use = "wilcox",
                          logfc.threshold = 0)

# 8. Add gene names as a column and filter for significant results
de_results$gene <- rownames(de_results)
de_results$neg_log10_pval <- -log10(de_results$p_val_adj)

# Filter for significant genes (adjust thresholds as needed)
sig_genes <- de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1, ]

# 9. Create volcano plot
ggplot(de_results, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1, ],
             color = "red", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linestyle = "dashed") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differential Expression: 24m vs 3m (Top 50% Prox Score)") +
  theme_minimal()

# 10. Extract significantly upregulated and downregulated genes
CB_assoc_genes <- sig_genes[sig_genes$avg_log2FC > 1, ]$gene  # Upregulated in 24m
CB_far_genes <- sig_genes[sig_genes$avg_log2FC < -1, ]$gene   # Downregulated in 24m

# Print results
cat("Upregulated genes in 24m vs 3m:\n")
print(head(CB_assoc_genes, 10))
cat("\nDownregulated genes in 24m vs 3m:\n") 
print(head(CB_far_genes, 10))

# View the DE results
head(de_results[order(de_results$p_val_adj), ], 10)

top_sig_genes <- head(de_results[order(de_results$p_val_adj), ], 10)$gene

# Create a combined plot for top genes
top_genes <- head(de_results[order(de_results$p_val_adj), ], 12)$gene
top_genes <- top_genes[top_genes %in% rownames(CB_gran_assoc_filtered)]

if(length(top_genes) > 0) {
  p_combined <- VlnPlot(CB_gran_assoc_filtered, 
                        features = top_genes, 
                        group.by = "Age",
                        pt.size = 0.1,
                        ncol = 4,
                        cols = c("3m" = "lightblue", "24m" = "coral"))
  
  ggsave("top_12_genes_violin_combined.pdf", p_combined, width = 24, height = 16)
  print(p_combined)
}


#12. Define gene lists and create module scores for CB_gran_assoc
LDAM <- c('Slc25a5', 'Npl', 'Angptl7','Pde2a','Ldhb','Cd63','Sepp1','Sdcbp','Adipor1','Rbbp4','Cndp2','Hsd17b4','Gpd11','Dazap2','Hnmpk','Rapsn','Cat','Kl','Nampt','Acsl1','Dpyd','Cd163')
neurop_full <- c('Axl','Cd9','Csf1r','Hif1a','Itgax','Tmem163','Apoe','Cybb','Lilr4b','Lgals3')

# Add module scores to CB_gran_assoc
CB_gran_assoc <- AddModuleScore(CB_gran_assoc, features = list(LDAM), name = "LDAM")
CB_gran_assoc <- AddModuleScore(CB_gran_assoc, features = list(neurop_full), name = "resilience")

# Note: AddModuleScore adds "1" suffix, so scores will be "LDAM1" and "resilience1"

# Function to create violin plots with statistics and strip plots
violin_with_stats <- function(seurat_obj, score_col, group_col, title_prefix, save_prefix) {
  library(ggplot2)
  library(dplyr)
  
  # Extract data
  plot_data <- seurat_obj@meta.data %>%
    select(all_of(c(score_col, group_col))) %>%
    na.omit()
  
  # Get unique groups
  groups <- unique(plot_data[[group_col]])
  
  # 1. Create violin plot
  p_violin <- VlnPlot(seurat_obj, 
                      features = score_col, 
                      group.by = group_col,
                      pt.size = 0.1,
                      cols = c("3m" = "lightblue", "24m" = "coral")) +
    ggtitle(paste(title_prefix, "by Age")) +
    xlab("Age") + ylab("Module Score") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(save_prefix, "_violin_plot.pdf"), p_violin, width = 8, height = 6)
  print(p_violin)
  
  # 2. Perform pairwise statistics (Wilcoxon test)
  if(length(groups) == 2) {
    group1_data <- plot_data[plot_data[[group_col]] == groups[1], score_col]
    group2_data <- plot_data[plot_data[[group_col]] == groups[2], score_col]
    
    wilcox_result <- wilcox.test(group1_data, group2_data)
    
    cat("\nPairwise Wilcoxon test results for", title_prefix, ":\n")
    cat(paste(groups[1], "vs", groups[2], ": p-value =", format(wilcox_result$p.value, scientific = TRUE)), "\n")
  }
  
  # 3. Create strip plot with regression line
  # Convert group to numeric for regression
  plot_data$group_numeric <- as.numeric(as.factor(plot_data[[group_col]]))
  
  # Fit linear model
  lm_fit <- lm(plot_data[[score_col]] ~ plot_data$group_numeric)
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  p_strip <- ggplot(plot_data, aes_string(x = group_col, y = score_col)) +
    geom_jitter(width = 0.2, alpha = 0.5, color = "gray") +
    geom_smooth(aes(x = group_numeric, y = get(score_col)), 
                method = "lm", se = TRUE, color = "red") +
    labs(title = paste("Strip Plot with Best Fit Line and 95% CI -", title_prefix),
         x = "Age", y = "Module Score") +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),
             hjust = 1.1, vjust = 1.1, color = "red") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(save_prefix, "_strip_with_line.pdf"), p_strip, width = 8, height = 6)
  print(p_strip)
  
  return(list(violin_plot = p_violin, strip_plot = p_strip, wilcox_result = wilcox_result))
}

CB_gran_assoc_filtered <- subset(CB_gran_assoc, Age %in% c("3", "24"))

# Create plots for LDAM score
violin_with_stats(CB_gran_assoc_filtered, "LDAM1", "Age", "LDAM Score", "LDAM")

# Create plots for resilience score  
violin_with_stats(CB_gran_assoc_filtered, "resilience1", "Age", "Resilience Score", "resilience")

