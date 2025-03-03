---
title: "Microglia Gene Module Analysis"
author: "James Haberberger"
---

# hdWGCNA Setup

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
theme_set(theme_cowplot())

set.seed(12345)
allowWGCNAThreads(nThreads = 5)

seurat_obj <- readRDS("data/microglia-aging-21MAR2023.rds")

seurat_obj <-
  NormalizeData(seurat_obj,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

seurat_obj <-
  FindVariableFeatures(seurat_obj,
                       selection.method = "vst",
                       nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <-
  RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Microglia" 
)

seurat_obj <- SetDatExpr(
  seurat_obj,
  use_metacells = F,
  assay="RNA"
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

plot_list <- PlotSoftPowers(seurat_obj)
power_table <- GetPowerTable(seurat_obj)

wrap_plots(plot_list, ncol=2)

# Network Construction

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Microglia' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='Microglia')

# Grab the Module Eigengenes, harmonize, plot major contributors.

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj)
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
seurat_obj <- ModuleConnectivity(seurat_obj, group_name = 'Microglia')
seurat_obj <- ResetModuleNames(seurat_obj, new_name = "Microglia-M")
PlotKMEs(seurat_obj, ncol=5)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

# Perform Module - Trait Correlation Analyses

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(region_coded = case_when(
    region == "CB" ~ 0,
    region == "TH" ~ 1,
    region == "HTH" ~ 2,
    region == "STR" ~ 3,
    region == "HP" ~ 4,
    region == "CTX" ~ 5,
  ))
# Gives an Error that ages were converted from factors in order 01, 03, and 22
# We are happy with that. The region has been coded in an order that makes sense
# given the dataset and our expectation of how disease progresses. 
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj, 
  traits = c("age", "region_coded"),
  features = "hMEs"
  )

p <- PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'orange',
  low_color = 'yellow',
  plot_max = 0.2,
  combine=TRUE,
  wgcna_name = "Microglia"
)
save_plot("output/module_trait_correlation.pdf", p, base_width = 11, base_height = 8.5)

# Go Term Enrichment Analysis

dbs <-
  c(
    'GO_Biological_Process_2021',
    'GO_Cellular_Component_2021',
    'GO_Molecular_Function_2021'
  )

seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs = dbs,
  max_genes = Inf,
  wait = T,
  wait_time = 1
)

enrich_df <- GetEnrichrTable(seurat_obj)

p <- EnrichrDotPlot(
  seurat_obj,
  mods = "all",
  database = "GO_Biological_Process_2021",
  n_terms = 1
)

save_plot(
  "output/gene_module_processes.pdf",
  plot = p,
  base_height = 8.5,
  base_width = 11
)
write.csv(enrich_df, "gene_module_processes.csv")

# Correlogram

ModuleCorrelogram(seurat_obj)

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)

mods <- levels(modules$module); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, features=mods, group.by = 'region')
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

p

plot_list = list()
model_output = list()

for (x in names(MEs)) {
  model_output[[x]] <- lme4::lmer(
    as.formula(paste0("`", x, "`", " ~ age + region + (1 + region|age) + (1 + age|region)")), 
    data=seurat_obj@meta.data
  ) 
  plot_list[[x]] <-  sjPlot::plot_model(
    model_output[[x]],
    show.values = TRUE,
    show.p = TRUE,
    sort.est = TRUE
  )
}

for (x in names(model_output)) {
    file_name = paste("output/mixed_effects_model/", x, ".pdf", sep="")
    pdf(file_name)
    print(plot_list[[x]])
    dev.off()
    
    file_name = paste("output/mixed_effects_model/", x, "_table.html", sep="")
    print(sjPlot::tab_model(model_output[[x]], file = file_name))
}

library(variancePartition)

vp <- fitExtractVarPartModel(
  seurat_obj@meta.data %>% 
    select(contains("Microglia-")) %>% 
    t(),
  ~ (1 | age) + (1 | region),
  seurat_obj@meta.data %>% 
    select(region, age) %>% 
    mutate(Individual = rownames(seurat_obj@meta.data))
)

vp <-
  vp[c(
    'Microglia-M1',
    'Microglia-M2',
    'Microglia-M3',
    'Microglia-M4',
    'Microglia-M5',
    'Microglia-M6',
    'Microglia-M7',
    'Microglia-M8',
    'Microglia-M9',
    'Microglia-M10',
    'Microglia-M11',
    'Microglia-M12',
    'Microglia-M13',
    'Microglia-M14',
    'Microglia-M15'
  ),]
p <- plotPercentBars(vp)
save_plot("output/variance_partition_barplot.pdf", p)
p <- plotVarPart(vp)
save_plot("output/variance_partition_violinplot.pdf", p)

# Get Region Based Expression

geneExpr <-
  seurat_obj@meta.data %>% select(contains("Microglia-")) %>% select(
    c(
      'Microglia-M1',
      'Microglia-M2',
      'Microglia-M3',
      'Microglia-M4',
      'Microglia-M5',
      'Microglia-M6',
      'Microglia-M7',
      'Microglia-M8',
      'Microglia-M9',
      'Microglia-M10',
      'Microglia-M11',
      'Microglia-M12',
      'Microglia-M13',
      'Microglia-M14',
      'Microglia-M15'
    )
  ) %>% t()




for (i in 1:15) {
  GE <-
    data.frame(
      Expression = geneExpr[i, ],
      Tissue = seurat_obj@meta.data %>% select(region, age) %>% mutate(Individual = rownames(seurat_obj@meta.data)) %>% pull(region)
    )
  p <-
    plotStratify(Expression ~ Tissue, GE, main = rownames(geneExpr)[i])
  save_plot(
    paste0(
      "output/gene_module_stratification/",
      i,
      "_region_stratification.pdf"
    ),
    p
  )
}

for (i in 1:15) {
  GE <-
    data.frame(
      Expression = geneExpr[i, ],
      Tissue = seurat_obj@meta.data %>% select(region, age) %>% mutate(Individual = rownames(seurat_obj@meta.data)) %>% pull(region),
      Age = seurat_obj@meta.data %>% select(region, age) %>% mutate(Individual = rownames(seurat_obj@meta.data)) %>% pull(age)
    )
  
  p <-
    plotStratify(Expression ~ Age, GE, main = rownames(geneExpr)[i])
  save_plot(
    paste0(
      "output/gene_module_stratification/",
      i,
      "_age_stratification.pdf"
    ),
    p
  )
}

for (i in 1:15) {
  GE <-
    data.frame(
      Expression = geneExpr[i, ],
      Tissue = seurat_obj@meta.data %>% select(region, age) %>% mutate(Individual = rownames(seurat_obj@meta.data)) %>% pull(region),
      Age = seurat_obj@meta.data %>% select(region, age) %>% mutate(Individual = rownames(seurat_obj@meta.data)) %>% pull(age)
    )
  
  p <-
    plotStratify(Expression ~ Age, GE, main = rownames(geneExpr)[i])
  save_plot(
    paste0(
      "output/gene_module_stratification/",
      i,
      "_age_stratification.pdf"
    ),
    p
  )
}