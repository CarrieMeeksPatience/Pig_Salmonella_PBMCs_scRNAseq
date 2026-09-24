{library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(msigdbr)
library(org.Hs.eg.db)
library(RColorBrewer)}

Monocyte_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_all <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")

Monocyte_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")

#Monocyte_GO_terms <- c("regulation of response to stimulus","immune system process","immune response","endocytosis","cell migration","cytokine production","interferon-mediated signaling pathway","positive regulation of cell activation","positive regulation of defense response")
Monocyte_GO_terms <- c("regulation of response to stimulus","positive regulation of defense response","positive regulation of cell activation","interferon−mediated signaling pathway","immune system process","immune response","endocytosis","cytokine production","cell migration")

# Ensure Description columns are character type
Monocyte_GO_terms<- trimws(as.character(Monocyte_GO_terms))
Monocyte_all$Description <- trimws(as.character(Monocyte_all$Description))
Monocyte_simp$Description <- trimws(as.character(Monocyte_simp$Description))
Monocyte_simp$Comparison <- gsub("D2 vs D0", "2 DPI vs 0 DPI", Monocyte_simp$Comparison)
Monocyte_simp$Comparison <- gsub("D8 vs D2", "8 DPI vs 2 DPI", Monocyte_simp$Comparison)
Monocyte_simp$Comparison <- gsub("D8 vs D0", "8 DPI vs 0 DPI", Monocyte_simp$Comparison)

Monocyte_simp_filt <- Monocyte_simp %>% filter(Description %in% Monocyte_GO_terms)

Monocyte_simp_filt <- Monocyte_simp_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Comparisons <- unique(Monocyte_simp_filt$Comparison)
# Set the desired order for the Comparison factor
Monocyte_simp_filt$Comparison <- factor(Monocyte_simp_filt$Comparison, levels = c("2 DPI vs 0 DPI", "8 DPI vs 2 DPI", "8 DPI vs 0 DPI"))
Monocyte_simp_filt_min <- min(Monocyte_simp_filt$Count)
Monocyte_simp_filt_max <- max(Monocyte_simp_filt$Count)
# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_Plots/Monocytes_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_04_23.pdf", width = 10, height = 10)
ggplot(Monocyte_simp_filt, aes(x = Comparison, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(Monocyte_simp_filt_min, Monocyte_simp_filt_max, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black", angle = 315, hjust =0),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Comparison", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in Monocytes (for paper)") 
dev.off()



NKcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_all <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")

NKcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")

NKcell_GO_terms <- c("response to stimulus","response to other organism","response to biotic stimulus","leukocyte activation","immune system process","immune response", "cell migration","apoptotic process")
# Ensure Description columns are character type
NKcell_GO_terms<- trimws(as.character(NKcell_GO_terms))
NKcell_all$Description <- trimws(as.character(NKcell_all$Description))
NKcell_simp$Description <- trimws(as.character(NKcell_simp$Description))
NKcell_simp$Comparison <- gsub("D2 vs D0", "2 DPI vs 0 DPI", NKcell_simp$Comparison)
NKcell_simp$Comparison <- gsub("D8 vs D2", "8 DPI vs 2 DPI", NKcell_simp$Comparison)
NKcell_simp$Comparison <- gsub("D8 vs D0", "8 DPI vs 0 DPI", NKcell_simp$Comparison)

NKcell_all$Comparison <- gsub("D2 vs D0", "2 DPI vs 0 DPI", NKcell_all$Comparison)
NKcell_all$Comparison <- gsub("D8 vs D2", "8 DPI vs 2 DPI", NKcell_all$Comparison)
NKcell_all$Comparison <- gsub("D8 vs D0", "8 DPI vs 0 DPI", NKcell_all$Comparison)


NKcell_simp_filt <- NKcell_simp %>% filter(Description %in% NKcell_GO_terms)

NKcell_simp_filt <- NKcell_simp_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Comparisons <- unique(NKcell_simp_filt$Comparison)
# Set the desired order for the Comparison factor
NKcell_simp_filt$Comparison <- factor(NKcell_simp_filt$Comparison, levels = c("2 DPI vs 0 DPI", "8 DPI vs 2 DPI", "8 DPI vs 0 DPI"))
NKcell_simp_filt_min <- min(NKcell_simp_filt$Count)
NKcell_simp_filt_max <- max(NKcell_simp_filt$Count)

# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_Plots/NKcells_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_04_23.pdf", width = 8, height = 8)
ggplot(NKcell_simp_filt, aes(x = Comparison, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(NKcell_simp_filt_min, NKcell_simp_filt_max, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black", angle = 315, hjust =0),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Comparison", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in NK Cells (for paper)") 
dev.off()

NKcell_all_filt <- NKcell_all %>% filter(Description %in% NKcell_GO_terms)

NKcell_all_filt <- NKcell_all_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Comparisons <- unique(NKcell_all_filt$Comparison)
# Set the desired order for the Comparison factor
NKcell_all_filt$Comparison <- factor(NKcell_all_filt$Comparison, levels = c("2 DPI vs 0 DPI", "8 DPI vs 2 DPI", "8 DPI vs 0 DPI"))
NKcell_all_filt_min <- min(NKcell_all_filt$Count)
NKcell_all_filt_max <- max(NKcell_all_filt$Count)
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_Plots/NKcells_pairwise_GOall_pos_BPinteresting_analysis_dotplot_short_2026_04_23.pdf",  width = 10, height = 8)
ggplot(NKcell_all_filt, aes(x = Comparison, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(NKcell_all_filt_min, NKcell_all_filt_max, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black", angle = 315, hjust =0),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Comparison", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in NK Cells (for paper)") 
dev.off()


Bcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_all <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")

Bcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")

Bcell_GO_terms <- c("response to other organism","positive regulation of response to stimulus","positive regulation of immune response","leukocyte activation","intracellular signal transduction","immune system process","immune response","apoptotic process")
# Ensure Description columns are character type
Bcell_GO_terms<- trimws(as.character(Bcell_GO_terms))
Bcell_all$Description <- trimws(as.character(Bcell_all$Description))
Bcell_simp$Description <- trimws(as.character(Bcell_simp$Description))
Bcell_simp$Comparison <- gsub("D2 vs D0", "2 DPI vs 0 DPI", Bcell_simp$Comparison)
Bcell_simp$Comparison <- gsub("D8 vs D2", "8 DPI vs 2 DPI", Bcell_simp$Comparison)
Bcell_simp$Comparison <- gsub("D8 vs D0", "8 DPI vs 0 DPI", Bcell_simp$Comparison)

Bcell_simp <- Bcell_simp %>% filter(Description %in% Bcell_GO_terms)
table(Bcell_simp$Description)
# Ensure Description columns are character type
Bcell_GO_terms<- trimws(as.character(Bcell_GO_terms))
Bcell_all$Description <- trimws(as.character(Bcell_all$Description))
Bcell_simp$Description <- trimws(as.character(Bcell_simp$Description))

Bcell_simp_filt <- Bcell_simp %>% filter(Description %in% Bcell_GO_terms)

Bcell_simp_filt <- Bcell_simp_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Comparisons <- unique(Bcell_simp_filt$Comparison)
# Set the desired order for the Comparison factor
Bcell_simp_filt$Comparison <- factor(Bcell_simp_filt$Comparison, levels = c("2 DPI vs 0 DPI", "8 DPI vs 2 DPI", "8 DPI vs 0 DPI"))
Bcell_simp_filt_min <- min(Bcell_simp_filt$Count)
Bcell_simp_filt_max <- max(Bcell_simp_filt$Count)
# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_Plots/Bcells_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_04_23.pdf", width = 10, height = 10)
ggplot(Bcell_simp_filt, aes(x = Comparison, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(Bcell_simp_filt_min, Bcell_simp_filt_max, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black", angle = 315, hjust =0),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Comparison", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in B Cells (for paper)") 
dev.off()



gdTcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_all <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")

gdTcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")

gdTcell_GO_terms <- c("actin binding","actin cytoskeleton organization","apoptotic process","cell adhesion molecule binding","cellular response to endogenous stimulus","leukocyte chemotaxis","lysosome","lytic vacuole","regulation of cell development","supramolecular fiber organization")
# Ensure Description columns are character type
gdTcell_GO_terms<- trimws(as.character(gdTcell_GO_terms))
gdTcell_all$Description <- trimws(as.character(gdTcell_all$Description))
gdTcell_simp$Description <- trimws(as.character(gdTcell_simp$Description))
gdTcell_all$Comparison <- gsub("D2 vs D0", "2 DPI vs 0 DPI", gdTcell_all$Comparison)
gdTcell_all$Comparison <- gsub("D8 vs D2", "8 DPI vs 2 DPI", gdTcell_all$Comparison)
gdTcell_all$Comparison <- gsub("D8 vs D0", "8 DPI vs 0 DPI", gdTcell_all$Comparison)


gdTcell_all_filt <- gdTcell_all %>% filter(Description %in% gdTcell_GO_terms)

gdTcell_all_filt <- gdTcell_all_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Comparisons <- unique(gdTcell_all_filt$Comparison)
# Set the desired order for the Comparison factor
gdTcell_all_filt$Comparison <- factor(gdTcell_all_filt$Comparison, levels = c("2 DPI vs 0 DPI", "8 DPI vs 2 DPI", "8 DPI vs 0 DPI"))
gdTcell_all_filt_min <- min(gdTcell_all_filt$Count)
gdTcell_all_filt_max <- max(gdTcell_all_filt$Count)
# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_Plots/CD2neg_gdTcell_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_04_23.pdf", width = 10, height = 12)
ggplot(gdTcell_all_filt, aes(x = Comparison, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(gdTcell_all_filt_min, gdTcell_all_filt_max, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black", angle = 315, hjust =0),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Comparison", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in CD2- gd T cells (for paper)") 
dev.off()


sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux 9.6 (Plow)

# Matrix products: default
# BLAS/LAPACK: /micromamba/envs/seurat+milo/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# time zone: America/Chicago
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#  [1] RColorBrewer_1.1-3     org.Hs.eg.db_3.18.0    AnnotationDbi_1.64.1  
#  [4] IRanges_2.36.0         S4Vectors_0.40.2       Biobase_2.62.0        
#  [7] BiocGenerics_0.48.1    msigdbr_2023.1.1       lubridate_1.9.4       
# [10] forcats_1.0.0          stringr_1.5.1          purrr_1.0.4           
# [13] readr_2.1.5            tidyr_1.3.1            tibble_3.2.1          
# [16] tidyverse_2.0.0        enrichplot_1.22.0      clusterProfiler_4.17.0
# [19] scCustomize_3.0.1      dplyr_1.1.4            ggplot2_3.5.2         
# [22] Matrix_1.6-5           Seurat_5.2.1           SeuratObject_5.0.2    
# [25] sp_2.2-0              

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.3.3           later_1.4.2            
#   [4] ggplotify_0.1.2         bitops_1.0-9            polyclip_1.10-7        
#   [7] janitor_2.2.1           fastDummies_1.7.5       lifecycle_1.0.4        
#  [10] globals_0.16.3          lattice_0.22-6          MASS_7.3-60            
#  [13] magrittr_2.0.3          plotly_4.10.4           httpuv_1.6.15          
#  [16] sctransform_0.4.1       spam_2.11-1             spatstat.sparse_3.1-0  
#  [19] reticulate_1.41.0.1     cowplot_1.1.3           pbapply_1.7-2          
#  [22] DBI_1.2.3               abind_1.4-8             zlibbioc_1.48.2        
#  [25] Rtsne_0.17              ggraph_2.2.1            RCurl_1.98-1.17        
#  [28] yulab.utils_0.2.0       tweenr_2.0.3            circlize_0.4.16        
#  [31] GenomeInfoDbData_1.2.11 ggrepel_0.9.6           irlba_2.3.5.1          
#  [34] listenv_0.9.1           spatstat.utils_3.1-2    tidytree_0.4.6         
#  [37] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-2  
#  [40] fitdistrplus_1.2-2      parallelly_1.42.0       codetools_0.2-20       
#  [43] ggforce_0.4.2           DOSE_3.28.2             tidyselect_1.2.1       
#  [46] shape_1.4.6.1           aplot_0.2.5             farver_2.1.2           
#  [49] viridis_0.6.5           matrixStats_1.5.0       spatstat.explore_3.3-4 
#  [52] jsonlite_2.0.0          tidygraph_1.3.1         progressr_0.15.1       
#  [55] ggridges_0.5.6          survival_3.8-3          tools_4.3.3            
#  [58] treeio_1.26.0           ica_1.0-3               Rcpp_1.0.14            
#  [61] glue_1.8.0              gridExtra_2.3           qvalue_2.34.0          
#  [64] GenomeInfoDb_1.38.8     withr_3.0.2             fastmap_1.2.0          
#  [67] digest_0.6.37           gridGraphics_0.5-1      timechange_0.3.0       
#  [70] R6_2.6.1                mime_0.13               ggprism_1.0.5          
#  [73] colorspace_2.1-1        scattermore_1.2         GO.db_3.18.0           
#  [76] tensor_1.5              spatstat.data_3.1-4     RSQLite_2.3.9          
#  [79] generics_0.1.3          data.table_1.17.0       graphlayouts_1.2.2     
#  [82] httr_1.4.7              htmlwidgets_1.6.4       scatterpie_0.2.4       
#  [85] uwot_0.2.3              pkgconfig_2.0.3         gtable_0.3.6           
#  [88] blob_1.2.4              lmtest_0.9-40           XVector_0.42.0         
#  [91] shadowtext_0.1.4        htmltools_0.5.8.1       dotCall64_1.2          
#  [94] fgsea_1.28.0            scales_1.4.0            png_0.1-8              
#  [97] spatstat.univar_3.1-2   snakecase_0.11.1        ggfun_0.1.8            
# [100] tzdb_0.4.0              reshape2_1.4.4          nlme_3.1-167           
# [103] zoo_1.8-13              cachem_1.1.0            GlobalOptions_0.1.2    
# [106] KernSmooth_2.23-26      parallel_4.3.3          miniUI_0.1.1.1         
# [109] vipor_0.4.7             HDO.db_0.99.1           ggrastr_1.0.2          
# [112] pillar_1.10.2           grid_4.3.3              vctrs_0.6.5            
# [115] RANN_2.6.2              promises_1.3.2          xtable_1.8-4           
# [118] cluster_2.1.8.1         beeswarm_0.4.0          paletteer_1.6.0        
# [121] cli_3.6.5               compiler_4.3.3          rlang_1.1.6            
# [124] crayon_1.5.3            future.apply_1.11.3     labeling_0.4.3         
# [127] rematch2_2.1.2          plyr_1.8.9              fs_1.6.6               
# [130] ggbeeswarm_0.7.2        stringi_1.8.7           viridisLite_0.4.2      
# [133] deldir_2.0-4            BiocParallel_1.36.0     babelgene_22.9         
# [136] Biostrings_2.70.3       lazyeval_0.2.2          spatstat.geom_3.3-5    
# [139] GOSemSim_2.28.1         RcppHNSW_0.6.0          hms_1.1.3              
# [142] patchwork_1.3.0         bit64_4.6.0-1           future_1.34.0          
# [145] KEGGREST_1.42.0         shiny_1.10.0            ROCR_1.0-11            
# [148] igraph_2.1.4            memoise_2.0.1           ggtree_3.10.1          
# [151] fastmatch_1.1-6         bit_4.6.0               gson_0.1.0             
# [154] ape_5.8-1              