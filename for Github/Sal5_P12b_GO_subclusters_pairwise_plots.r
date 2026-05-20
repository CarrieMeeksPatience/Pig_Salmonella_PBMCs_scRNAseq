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


Monocyte_all <-read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_Monocytes_subcluster_pairwise_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
#rename "Comparison" to "Cluster"
colnames(Monocyte_all)[12] <- "Cluster"
# use Gsub to replace "0_1 vs 0_0" with "MC1"
Monocyte_all$Cluster <- gsub("0_1 vs 0_0", "MC1", Monocyte_all$Cluster)
Monocyte_all$Cluster <- gsub("0_0 vs 0_1", "MC0", Monocyte_all$Cluster)

Monocyte_GO_terms_old <- c("response to molecule of bacterial origin","response to lipopolysaccharide","positive regulation of inflammatory response","positive regulation of apoptotic process","myeloid cell differentiation","innate immune response−activating signaling pathway","endocytosis","cytokine production","cellular response to type II interferon","cellular response to type I interferon","cellular response to lipopolysaccharide","cell development","activation of immune response")

# Ensure Description columns are character type
Monocyte_GO_terms_old <- trimws(as.character(Monocyte_GO_terms_old))
Monocyte_all$Description <- trimws(as.character(Monocyte_all$Description))

Monocyte_all_filt <- Monocyte_all %>% filter(Description %in% Monocyte_GO_terms_old)

Monocyte_all_filt <- Monocyte_all_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Clusters <- unique(Monocyte_all_filt$Cluster)
# Set the desired order for the Comparison factor
Monocyte_all_filt$Cluster <- factor(Monocyte_all_filt$Cluster, levels = c("MC0", "MC1"))
# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_Plots/Monocytes_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_05_05.pdf", width = 10, height = 10)
ggplot(Monocyte_all_filt, aes(x = Cluster, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(9, 67, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black",  hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Cluster", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in Monocytes (for paper)") 
dev.off()


NKcells_all <-read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_NKcells_subcluster_pairwise_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
#rename "Comparison" to "Cluster"
colnames(NKcells_all)[12] <- "Cluster"
# use Gsub to replace "0_1 vs 0_0" with "NKcell 1"
NKcells_all$Cluster <- gsub("4_1 vs 4_0", "NKC1", NKcells_all$Cluster)
NKcells_all$Cluster <- gsub("4_0 vs 4_1", "NKC0", NKcells_all$Cluster)

NKcells_GO_terms_old <- c("regulation of cytokine production","positive regulation of leukocyte mediated cytotoxicity","positive regulation of cell migration","positive regulation of cell killing","natural killer cell mediated cytotoxicity","natural killer cell activation","leukocyte mediated immunity","leukocyte mediated cytotoxicity","leukocyte activation involved in immune response","immune response","immune effector process","cytokine−mediated signaling pathway","cellular response to cytokine stimulus","cell population proliferation")

# Ensure Description columns are character type
NKcells_GO_terms_old <- trimws(as.character(NKcells_GO_terms_old))
NKcells_all$Description <- trimws(as.character(NKcells_all$Description))

NKcells_all_filt <- NKcells_all %>% filter(Description %in% NKcells_GO_terms_old)

NKcells_all_filt <- NKcells_all_filt %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Clusters <- unique(NKcells_all_filt$Cluster)
# Set the desired order for the Comparison factor
NKcells_all_filt$Cluster <- factor(NKcells_all_filt$Cluster, levels = c("NKC0", "NKC1"))

# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_Plots/NKcells_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_05_05.pdf", width = 10, height = 10)
ggplot(NKcells_all_filt, aes(x = Cluster, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(20, 71, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black",  hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Cluster", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in NK cells (for paper)") 
dev.off()

NKcells_GO_terms_new <- c("actin cytoskeleton organization","cell development","positive regulation of T cell activation","leukocyte migration","leukocyte activation involved in immune response","immune response","immunological synapse","natural killer cell mediated cytotoxicity","leukocyte mediated cytotoxicity","positive regulation of cytokine production","killing of cells of another organism","cytokine production","actin filament organization","immune response")
# Ensure Description columns are character type
NKcells_GO_terms_new <- trimws(as.character(NKcells_GO_terms_new))
NKcells_all$Description <- trimws(as.character(NKcells_all$Description))

NKcells_all_filt2 <- NKcells_all %>% filter(Description %in% NKcells_GO_terms_new)

NKcells_all_filt2 <- NKcells_all_filt2 %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Clusters <- unique(NKcells_all_filt2$Cluster)
# Set the desired order for the Comparison factor
NKcells_all_filt2$Cluster <- factor(NKcells_all_filt2$Cluster, levels = c("NKC0", "NKC1"))

# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_Plots/NKcells_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_05_05_V2.pdf", width = 10, height = 12)
ggplot(NKcells_all_filt2, aes(x = Cluster, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(20, 71, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black",  hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Cluster", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in NK cells (for paper)") 
dev.off()

NKcells_GO_terms_new_short <- c("actin cytoskeleton organization","cell development","leukocyte migration","leukocyte activation involved in immune response","immune response","natural killer cell mediated cytotoxicity","leukocyte mediated cytotoxicity","positive regulation of cytokine production","killing of cells of another organism","immune response","lymphocyte proliferation","positive regulation of cell differentiation")
# Ensure Description columns are character type
NKcells_GO_terms_new_short<- trimws(as.character(NKcells_GO_terms_new_short))
NKcells_all$Description <- trimws(as.character(NKcells_all$Description))

NKcells_all_filt3 <- NKcells_all %>% filter(Description %in% NKcells_GO_terms_new_short)

NKcells_all_filt3 <- NKcells_all_filt3 %>% filter(p.adjust <= 0.05 & Count >= 5) # Filter to keep only GO terms with Count >= 20 & p.adjust <= 0.05
# Ensure Cluster is a factor and set levels to include all clusters
all_Clusters <- unique(NKcells_all_filt3$Cluster)
# Set the desired order for the Comparison factor
NKcells_all_filt3$Cluster <- factor(NKcells_all_filt3$Cluster, levels = c("NKC0", "NKC1"))
# Now make a dotplot, use GGplot2, split the plot by cluster
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_Plots/NKcells_pairwise_GOsimp_pos_BPinteresting_analysis_dotplot_short_2026_05_07_V3.pdf", width = 10, height = 11)
ggplot(NKcells_all_filt3, aes(x = Cluster, y = Description, size = Count, fill = p.adjust)) +
geom_point(shape = 21) +
    scale_fill_distiller(
    palette = "PuBuGn", 
    direction = 1, 
    trans = "reverse") +
scale_size(range = c(7, 28), breaks = seq(5, 80, length.out = 4)) + 
  theme_classic() + # Remove gray background
  theme(
   axis.text.x = element_text(size = 16, color = "black",  hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    panel.spacing = unit(0.01, "lines")  # Reduce space between columns
  ) +
  labs(x = "Cluster", y = "GO Term", fill = "Adjusted p-value", size = "Count") + 
  ggtitle("Interesting GO Terms from Pos DEGs in NK cells (for paper)") 
dev.off()


# Now do GO_DEG plots
Monocyte_1v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
# Get the number of significant positive and negative DEGs
pos_DEGs <- subset(Monocyte_1v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
neg_DEGs <- subset(Monocyte_1v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))

# Count the number of positive and negative DEGs
Monocyte_1v0_num_pos_DEGs <- nrow(pos_DEGs)
Monocyte_1v0_num_neg_DEGs <- nrow(neg_DEGs)
rm(pos_DEGs, neg_DEGs)
# Print the results
Monocyte_0v1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
# Get the number of significant positive and negative DEGs
pos_DEGs <- subset(Monocyte_0v1, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
neg_DEGs <- subset(Monocyte_0v1, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))

# Count the number of positive and negative DEGs
Monocyte_0v1_num_pos_DEGs <- nrow(pos_DEGs)
Monocyte_0v1_num_neg_DEGs <- nrow(neg_DEGs)
rm(pos_DEGs, neg_DEGs)
number_DEGs <- rbind(Monocyte_1v0_num_pos_DEGs, Monocyte_1v0_num_neg_DEGs, Monocyte_0v1_num_pos_DEGs, Monocyte_0v1_num_neg_DEGs)
write.table(number_DEGs, "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_res009_pairwise_number_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
Monocyte_DE_all <- rbind(Monocyte_1v0, Monocyte_0v1)
#rename "Comparison" to "Cluster"
colnames(Monocyte_DE_all)[9] <- "Cluster"
# use Gsub to replace "0_1 vs 0_0" with "MC1"
Monocyte_DE_all$Cluster <- gsub("0_1 vs 0_0", "MC1", Monocyte_DE_all$Cluster)
Monocyte_DE_all$Cluster <- gsub("0_0 vs 0_1", "MC0", Monocyte_DE_all$Cluster)
Monocyte_DE_all$pct.1 <- Monocyte_DE_all$pct.1 * 100
Monocyte_DE_all$pct.2 <- Monocyte_DE_all$pct.2 * 100
Monocyte_DE_all$"pct.1-pct.2"<- Monocyte_DE_all$pct.1 - Monocyte_DE_all$pct.2
ORG<- read.csv(file ="/Annotation_files/PigToHuman_GeneOrthos_v11_1_97_scGenes.csv", header = T,row.names=1)
colnames(ORG)[11] <- "Gene"
# Replace "_" with "-" in the Gene column
ORG$Gene <- gsub("_", "-", ORG$Gene)
# For duplicates in the Gene column, keep the row with the highest X.id..query.gene.identical.to.target.Human.gene
#ORG2 <- ORG %>% group_by(Gene) %>% filter(X.id..query.gene.identical.to.target.Human.gene == max(X.id..query.gene.identical.to.target.Human.gene)) %>% ungroup()
#subset the columns needed for annotation
ORG2.subset <- ORG %>% dplyr::select(c(1,11,7,8,10))
ORG2.subset$Gene <- as.character(ORG2.subset$Gene)
Monocyte_DE_all<-left_join(Monocyte_DE_all,ORG2.subset,by=c("gene"="Gene"))
Monocyte_DE_all <- unique(Monocyte_DE_all)
# Set the desired order for the Comparison factor
all_Comparisons <- unique(Monocyte_DE_all$Cluster)

#Make the old marker plot for the paper
Monocytes_markers_old <- c("HDAC9","ERBIN","AKT3","ANKRD17","NEK7","BIRC3","NFKBIA","IKBKB","DOCK4","ITGA4","APP","ITGAL","MERTK","CD74","S100A10","FCGR3A","TNFRSF1B","RACK1","BST2","B2M","CXCL10","GBP2","TGFB1","GRN","MYD88","STAT1","GBP1","IRF8","LYZ","STAT3","LTF","SOD2","CD163","MAPK14","S100A8","CD14","S100A9","FOS")

Monocyte_GO_DEGs_old <- Monocyte_DE_all %>%
  dplyr::filter(gene %in% Monocytes_markers_old)
  

write.table(Monocyte_GO_DEGs_old, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/Monocyte_GO_DEGs_old.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GO_DEGs_old_DEGs<-Monocyte_GO_DEGs_old$gene
Monocyte_GO_DEGs_old_DEGs <- unique(Monocyte_GO_DEGs_old_DEGs)
print("length(Monocyte_GO_DEGs_old_DEGs)")
length(Monocyte_GO_DEGs_old_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GO_DEGs_old$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GO_DEGs_old$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GO_DEGs_old_breaks <- c(
    (min(Monocyte_GO_DEGs_old$avg_log2FC) + (min(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GO_DEGs_old$avg_log2FC) + (max(Monocyte_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GO_DEGs_old_breaks<-round(Monocyte_GO_DEGs_old_breaks, 2)  # Round normally for the rest
Monocyte_GO_DEGs_old_breaks<- c(min_value, Monocyte_GO_DEGs_old_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GO_DEGs_old_breaks")
Monocyte_GO_DEGs_old_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
Monocyte_GO_DEGs_old_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
Monocyte_GO_DEGs_old_breaks<-round(Monocyte_GO_DEGs_old_breaks, 2)  # Round normally for the rest
print("Monocyte_GO_DEGs_old_breaks")
Monocyte_GO_DEGs_old_breaks
}
all_Clusters <- unique(Monocyte_GO_DEGs_old$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOmarkers_dotPlotCSP.pdf", width=6 ,height=14)
ggplot( Monocyte_GO_DEGs_old, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GO_DEGs_old_breaks,
      labels =  Monocyte_GO_DEGs_old_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in Monocyte Subclusters")  # Add title
dev.off()

#Make the new marker plot for the paper
Monocytes_markers_new <- c("HDAC9","ERBIN","AKT3","ANKRD17","NEK7","BIRC3","NFKBIA","IKBKB","DOCK4","ITGA4","APP","ITGAL","MERTK","CD74","S100A10","FCGR3A","TNFRSF1B","RACK1","BST2","B2M","CXCL10","GBP2","TGFB1","GRN","MYD88","STAT1","GBP1","IRF8","LYZ","STAT3","LTF","SOD2","CD163","MAPK14","S100A8","CD14","S100A9","FOS","CD52","SLA-DQA1","SLA-DQB1","SLA-DRA","SLA-DRB1","SIGLEC1","NFKB1")

Monocyte_GO_DEGs_new <- Monocyte_DE_all %>%
  dplyr::filter(gene %in% Monocytes_markers_new)
  

write.table(Monocyte_GO_DEGs_new, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/Monocyte_GO_DEGs_new.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GO_DEGs_new_DEGs<-Monocyte_GO_DEGs_new$gene
Monocyte_GO_DEGs_new_DEGs <- unique(Monocyte_GO_DEGs_new_DEGs)
print("length(Monocyte_GO_DEGs_new_DEGs)")
length(Monocyte_GO_DEGs_new_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GO_DEGs_new$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GO_DEGs_new$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GO_DEGs_new_breaks <- c(
    (min(Monocyte_GO_DEGs_new$avg_log2FC) + (min(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GO_DEGs_new$avg_log2FC) + (max(Monocyte_GO_DEGs_new$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GO_DEGs_new_breaks<-round(Monocyte_GO_DEGs_new_breaks, 2)  # Round normally for the rest
Monocyte_GO_DEGs_new_breaks<- c(min_value, Monocyte_GO_DEGs_new_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GO_DEGs_new_breaks")
Monocyte_GO_DEGs_new_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
Monocyte_GO_DEGs_new_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
Monocyte_GO_DEGs_new_breaks<-round(Monocyte_GO_DEGs_new_breaks, 2)  # Round normally for the rest
print("Monocyte_GO_DEGs_new_breaks")
Monocyte_GO_DEGs_new_breaks
}
all_Clusters <- unique(Monocyte_GO_DEGs_new$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOmarkers_dotPlotCSP.pdf", width=5 ,height=15)
ggplot( Monocyte_GO_DEGs_new, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GO_DEGs_new_breaks,
      labels =  Monocyte_GO_DEGs_new_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in Monocyte Subclusters")  # Add title
dev.off()

print("Monocyte response to lipopolysaccharide Plot")
response_to_lipopolysaccharide <- Monocyte_all %>% 
  dplyr::filter(Description == "response to lipopolysaccharide") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  Monocyte_GOrtl <- Monocyte_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% response_to_lipopolysaccharide)) %>%  # Check if any row satisfies the condition
      ungroup() 

write.table(Monocyte_GOrtl, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/Monocyte_GOrtl_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOrtl_DEGs<-Monocyte_GOrtl$gene
Monocyte_GOrtl_DEGs <- unique(Monocyte_GOrtl_DEGs)
print("length(Monocyte_GOrtl_DEGs)")
length(Monocyte_GOrtl_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOrtl$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOrtl$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOrtl_breaks <- c(
    (min(Monocyte_GOrtl$avg_log2FC) + (min(Monocyte_GOrtl$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOrtl$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOrtl$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOrtl$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOrtl$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOrtl$avg_log2FC) + (max(Monocyte_GOrtl$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOrtl_breaks<-round(Monocyte_GOrtl_breaks, 2)  # Round normally for the rest
Monocyte_GOrtl_breaks<- c(min_value, Monocyte_GOrtl_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOrtl_breaks")
Monocyte_GOrtl_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
Monocyte_GOrtl_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
Monocyte_GOrtl_breaks<-round(Monocyte_GOrtl_breaks, 2)  # Round normally for the rest
print("Monocyte_GOrtl_breaks")
Monocyte_GOrtl_breaks
}
all_Clusters <- unique(Monocyte_GOrtl$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOresponse_to_lipopolysaccharide_dotPlotCSP.pdf", width=5 ,height=11)
ggplot( Monocyte_GOrtl, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOrtl_breaks,
      labels =  Monocyte_GOrtl_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black",hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to lipopolysaccharide' in Monocytes")  # Add title
dev.off()


print("Monocyte myeloid cell differentiation Plot")
myeloid_cell_differentiation<- Monocyte_all %>% 
  dplyr::filter(Description == "myeloid cell differentiation") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
Monocyte_GOmcd <- Monocyte_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(pct.1 >= 50 | pct.2 >= 50, abs(avg_log2FC) >= 0.5), Human.gene.stable.ID %in% myeloid_cell_differentiation) %>%  # Check if any row satisfies the condition
      ungroup() 

#write.table(Monocyte_GOmcd, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/Monocyte_GOmcd_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOmcd_DEGs<-Monocyte_GOmcd$gene
Monocyte_GOmcd_DEGs <- unique(Monocyte_GOmcd_DEGs)
print("length(Monocyte_GOmcd_DEGs)")
length(Monocyte_GOmcd_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOmcd$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOmcd$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOmcd_breaks <- c(
    (min(Monocyte_GOmcd$avg_log2FC) + (min(Monocyte_GOmcd$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOmcd$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOmcd$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOmcd$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOmcd$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOmcd$avg_log2FC) + (max(Monocyte_GOmcd$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOmcd_breaks<-round(Monocyte_GOmcd_breaks, 2)  # Round normally for the rest
Monocyte_GOmcd_breaks<- c(min_value, Monocyte_GOmcd_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOmcd_breaks")
Monocyte_GOmcd_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
Monocyte_GOmcd_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
Monocyte_GOmcd_breaks<-round(Monocyte_GOmcd_breaks, 2)  # Round normally for the rest
print("Monocyte_GOmcd_breaks")
Monocyte_GOmcd_breaks
}
all_Clusters <- unique(Monocyte_GOmcd$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOmyeloid_cell_differentiation_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOmcd, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOmcd_breaks,
      labels =  Monocyte_GOmcd_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'myeloid cell differentiation' in Monocytes")  # Add title
dev.off()

print("Monocyte cytokine production Plot")
cytokine_production <- Monocyte_all %>% 
  dplyr::filter(Description == "cytokine production") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  Monocyte_GOcp <- Monocyte_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(pct.1 >= 50 | pct.2 >= 50, abs(avg_log2FC) >= 1),Human.gene.stable.ID %in% cytokine_production) %>%  # Check if any row satisfies the condition
      ungroup() 

#write.table(Monocyte_GOcp, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/Monocyte_GOcp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOcp_DEGs<-Monocyte_GOcp$gene
Monocyte_GOcp_DEGs <- unique(Monocyte_GOcp_DEGs)
print("length(Monocyte_GOcp_DEGs)")
length(Monocyte_GOcp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOcp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOcp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOcp_breaks <- c(
    (min(Monocyte_GOcp$avg_log2FC) + (min(Monocyte_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOcp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOcp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOcp$avg_log2FC) + (max(Monocyte_GOcp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOcp_breaks<-round(Monocyte_GOcp_breaks, 2)  # Round normally for the rest
Monocyte_GOcp_breaks<- c(min_value, Monocyte_GOcp_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOcp_breaks")
Monocyte_GOcp_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
Monocyte_GOcp_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
Monocyte_GOcp_breaks<-round(Monocyte_GOcp_breaks, 2)  # Round normally for the rest
print("Monocyte_GOcp_breaks")
Monocyte_GOcp_breaks
}
all_Clusters <- unique(Monocyte_GOcp$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOcytokine_production_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOcp, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOcp_breaks,
      labels =  Monocyte_GOcp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cytokine production' in Monocytes")  # Add title
dev.off()


NKcells_1v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
pos_DEGs <- subset(NKcells_1v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
neg_DEGs <- subset(NKcells_1v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))

# Count the number of positive and negative DEGs
NKcells_1v0_num_pos_DEGs <- nrow(pos_DEGs)
NKcells_1v0_num_neg_DEGs <- nrow(neg_DEGs)
rm(pos_DEGs, neg_DEGs)

NKcells_0v1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
# Get the number of significant positive and negative DEGs
pos_DEGs <- subset(NKcells_0v1, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
neg_DEGs <- subset(NKcells_0v1, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))

# Count the number of positive and negative DEGs
NKcells_0v1_num_pos_DEGs <- nrow(pos_DEGs)
NKcells_0v1_num_neg_DEGs <- nrow(neg_DEGs)
rm(pos_DEGs, neg_DEGs)
number_DEGs <- rbind(NKcells_1v0_num_pos_DEGs, NKcells_1v0_num_neg_DEGs, NKcells_0v1_num_pos_DEGs, NKcells_0v1_num_neg_DEGs)
write.table(number_DEGs, "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_res009_pairwise_number_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)

NKcells_DE_all <- rbind(NKcells_1v0, NKcells_0v1)
#rename "Comparison" to "Cluster"
colnames(NKcells_DE_all)[9] <- "Cluster"
# use Gsub to replace "4_1 vs 4_0" with "NKC1"
NKcells_DE_all$Cluster <- gsub("4_1 vs 4_0", "NKC1", NKcells_DE_all$Cluster)
NKcells_DE_all$Cluster <- gsub("4_0 vs 4_1", "NKC0", NKcells_DE_all$Cluster)
NKcells_DE_all$pct.1 <- NKcells_DE_all$pct.1 * 100
NKcells_DE_all$pct.2 <- NKcells_DE_all$pct.2 * 100
NKcells_DE_all$"pct.1-pct.2"<- NKcells_DE_all$pct.1 - NKcells_DE_all$pct.2
ORG<- read.csv(file ="/Annotation_files/PigToHuman_GeneOrthos_v11_1_97_scGenes.csv", header = T,row.names=1)
colnames(ORG)[11] <- "Gene"
# Replace "_" with "-" in the Gene column
ORG$Gene <- gsub("_", "-", ORG$Gene)
# For duplicates in the Gene column, keep the row with the highest X.id..query.gene.identical.to.target.Human.gene
#ORG2 <- ORG %>% group_by(Gene) %>% filter(X.id..query.gene.identical.to.target.Human.gene == max(X.id..query.gene.identical.to.target.Human.gene)) %>% ungroup()
#subset the columns needed for annotation
ORG2.subset <- ORG %>% dplyr::select(c(1,11,7,8,10))
ORG2.subset$Gene <- as.character(ORG2.subset$Gene)
NKcells_DE_all<-left_join(NKcells_DE_all,ORG2.subset,by=c("gene"="Gene"))
NKcells_DE_all <- unique(NKcells_DE_all)
# Set the desired order for the Comparison factor
all_Comparisons <- unique(NKcells_DE_all$Cluster)



print("NKcells natural killer cell mediated cytotoxicity Plot")
natural_killer_cell_mediated_cytotoxicity <- NKcells_all %>% 
  dplyr::filter(Description == "natural killer cell mediated cytotoxicity") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  NKcells_GOnkcmc <- NKcells_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% natural_killer_cell_mediated_cytotoxicity)) %>%  # Check if any row satisfies the condition
      ungroup() 

write.table(NKcells_GOnkcmc, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOnkcmc_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcells_GOnkcmc_DEGs<-NKcells_GOnkcmc$gene
NKcells_GOnkcmc_DEGs <- unique(NKcells_GOnkcmc_DEGs)
print("length(NKcells_GOnkcmc_DEGs)")
length(NKcells_GOnkcmc_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOnkcmc$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOnkcmc$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOnkcmc_breaks <- c(
    (min(NKcells_GOnkcmc$avg_log2FC) + (min(NKcells_GOnkcmc$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOnkcmc$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOnkcmc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOnkcmc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOnkcmc$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOnkcmc$avg_log2FC) + (max(NKcells_GOnkcmc$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOnkcmc_breaks<-round(NKcells_GOnkcmc_breaks, 2)  # Round normally for the rest
NKcells_GOnkcmc_breaks<- c(min_value, NKcells_GOnkcmc_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOnkcmc_breaks")
NKcells_GOnkcmc_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOnkcmc_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOnkcmc_breaks<-round(NKcells_GOnkcmc_breaks, 2)  # Round normally for the rest
print("NKcells_GOnkcmc_breaks")
NKcells_GOnkcmc_breaks
}
all_Clusters <- unique(NKcells_GOnkcmc$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOnatural killer_cell_mediated_cytotoxicity_dotPlotCSP.pdf", width=6 ,height=6)
ggplot( NKcells_GOnkcmc, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOnkcmc_breaks,
      labels =  NKcells_GOnkcmc_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0, size = 8, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'natural killer cell mediated cytotoxicity' in NK cells")  # Add title
dev.off()

print("NKcells leukocyte migration Plot")
leukocyte_migration <- NKcells_all %>% 
  dplyr::filter(Description == "leukocyte migration") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  NKcells_GOlm <- NKcells_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% leukocyte_migration)) %>%  # Check if any row satisfies the condition
      ungroup() 

write.table(NKcells_GOlm, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOlm_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcells_GOlm_DEGs<-NKcells_GOlm$gene
NKcells_GOlm_DEGs <- unique(NKcells_GOlm_DEGs)
print("length(NKcells_GOlm_DEGs)")
length(NKcells_GOlm_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOlm$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOlm$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOlm_breaks <- c(
    (min(NKcells_GOlm$avg_log2FC) + (min(NKcells_GOlm$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOlm$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOlm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOlm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOlm$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOlm$avg_log2FC) + (max(NKcells_GOlm$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOlm_breaks<-round(NKcells_GOlm_breaks, 2)  # Round normally for the rest
NKcells_GOlm_breaks<- c(min_value, NKcells_GOlm_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOlm_breaks")
NKcells_GOlm_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOlm_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOlm_breaks<-round(NKcells_GOlm_breaks, 2)  # Round normally for the rest
print("NKcells_GOlm_breaks")
NKcells_GOlm_breaks
}
all_Clusters <- unique(NKcells_GOlm$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOleukocyte_migration_dotPlotCSP.pdf", width=5 ,height=13)
ggplot( NKcells_GOlm, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOlm_breaks,
      labels =  NKcells_GOlm_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte migration' in NK cells")  # Add title
dev.off()

print("NKcells cell development Plot")
cell_development <- NKcells_all %>% 
  dplyr::filter(Description == "cell development") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  NKcells_GOcd <- NKcells_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(pct.1 >= 50 | pct.2 >= 50, abs(avg_log2FC) >= 1.5), Human.gene.stable.ID %in% cell_development) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcells_GOcd, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOcd_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
cell_development_paper <- c("CARD11","IKZF3","IL4R","ITPKB","PPP3CA","PTK2B","PTPRC","RORA","RUNX1","TGFBR2","TOX","WNK1","ZEB2")
  NKcells_GOcd_paper <- NKcells_DE_all %>%
dplyr::filter(gene %in% cell_development_paper) 

#make vector of genes
NKcells_GOcd_paper_DEGs<-NKcells_GOcd_paper$gene
NKcells_GOcd_paper_DEGs <- unique(NKcells_GOcd_paper_DEGs)
print("length(NKcells_GOcd_paper_DEGs)")
length(NKcells_GOcd_paper_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOcd_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOcd_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOcd_paper_breaks <- c(
    (min(NKcells_GOcd_paper$avg_log2FC) + (min(NKcells_GOcd_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOcd_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOcd_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOcd_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOcd_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOcd_paper$avg_log2FC) + (max(NKcells_GOcd_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOcd_paper_breaks<-round(NKcells_GOcd_paper_breaks, 2)  # Round normally for the rest
NKcells_GOcd_paper_breaks<- c(min_value, NKcells_GOcd_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOcd_paper_breaks")
NKcells_GOcd_paper_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOcd_paper_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOcd_paper_breaks<-round(NKcells_GOcd_paper_breaks, 2)  # Round normally for the rest
print("NKcells_GOcd_paper_breaks")
NKcells_GOcd_paper_breaks
}
all_Clusters <- unique(NKcells_GOcd_paper$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOcell_development_dotPlotCSP.pdf", width=5 ,height=7)
ggplot( NKcells_GOcd_paper, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOcd_paper_breaks,
      labels =  NKcells_GOcd_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cell development' in NK cells")  # Add title
dev.off()

print("NKcells cytokine production Plot")
cytokine_production <- NKcells_all %>% 
  dplyr::filter(Description == "cytokine production") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  NKcells_GOcp <- NKcells_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
 dplyr::filter(any(pct.1 >= 50 | pct.2 >= 50, abs(avg_log2FC) >= 0.75),Human.gene.stable.ID %in% cytokine_production) %>%  # Check if any row satisfies the condition
      ungroup() 

#write.table(NKcells_GOcp, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOcp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcells_GOcp_DEGs<-NKcells_GOcp$gene
NKcells_GOcp_DEGs <- unique(NKcells_GOcp_DEGs)
print("length(NKcells_GOcp_DEGs)")
length(NKcells_GOcp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOcp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOcp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOcp_breaks <- c(
    (min(NKcells_GOcp$avg_log2FC) + (min(NKcells_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOcp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOcp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOcp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOcp$avg_log2FC) + (max(NKcells_GOcp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOcp_breaks<-round(NKcells_GOcp_breaks, 2)  # Round normally for the rest
NKcells_GOcp_breaks<- c(min_value, NKcells_GOcp_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOcp_breaks")
NKcells_GOcp_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOcp_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOcp_breaks<-round(NKcells_GOcp_breaks, 2)  # Round normally for the rest
print("NKcells_GOcp_breaks")
NKcells_GOcp_breaks
}
all_Clusters <- unique(NKcells_GOcp$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOcytokine_production_dotPlotCSP.pdf", width=6 ,height=13)
ggplot( NKcells_GOcp, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOcp_breaks,
      labels =  NKcells_GOcp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cytokine production' in NK cells")  # Add title
dev.off()


#Make the old marker plot for the paper
NKcell_markers_old <- c("FCER1G","B2M","CD2","CD74","CORO1A","HCST","HSPA8","KLRK1","NKG7","RAC2","S100A13","TYROBP","ADGRG1","CARD11","CBLB","IKZF3","IL12RB2","ITPKB","RORA","RUNX1","TOX","ZEB2","PPP3CA","DOCK8","PTK2B","STAT3","SUN2","TGFBR2","WNK1","AKT3","PTPRC","ITGA4")

NKcell_GO_DEGs_old <- NKcells_DE_all %>%
  dplyr::filter(gene %in% NKcell_markers_old)
  

write.table(NKcell_GO_DEGs_old, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcell_GO_DEGs_old.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GO_DEGs_old_DEGs<-NKcell_GO_DEGs_old$gene
NKcell_GO_DEGs_old_DEGs <- unique(NKcell_GO_DEGs_old_DEGs)
print("length(NKcell_GO_DEGs_old_DEGs)")
length(NKcell_GO_DEGs_old_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GO_DEGs_old$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GO_DEGs_old$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GO_DEGs_old_breaks <- c(
    (min(NKcell_GO_DEGs_old$avg_log2FC) + (min(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GO_DEGs_old$avg_log2FC) + (max(NKcell_GO_DEGs_old$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GO_DEGs_old_breaks<-round(NKcell_GO_DEGs_old_breaks, 2)  # Round normally for the rest
NKcell_GO_DEGs_old_breaks<- c(min_value, NKcell_GO_DEGs_old_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GO_DEGs_old_breaks")
NKcell_GO_DEGs_old_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcell_GO_DEGs_old_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcell_GO_DEGs_old_breaks<-round(NKcell_GO_DEGs_old_breaks, 2)  # Round normally for the rest
print("NKcell_GO_DEGs_old_breaks")
NKcell_GO_DEGs_old_breaks
}
all_Clusters <- unique(NKcell_GO_DEGs_old$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOmarkers_dotPlotCSP.pdf", width=5 ,height=14)
ggplot( NKcell_GO_DEGs_old, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GO_DEGs_old_breaks,
      labels =  NKcell_GO_DEGs_old_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in NK cell Subclusters")  # Add title
dev.off()


print("NKcells leukocyte mediated cytotoxicity Plot")
leukocyte_mediated_cytotoxicity <- NKcells_all %>% 
  dplyr::filter(Description == "leukocyte mediated cytotoxicity") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  NKcells_GOlmc <- NKcells_DE_all %>%
 dplyr::filter(Human.gene.stable.ID %in% leukocyte_mediated_cytotoxicity) 

write.table(NKcells_GOlmc, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOlmc_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcells_GOlmc_DEGs<-NKcells_GOlmc$gene
NKcells_GOlmc_DEGs <- unique(NKcells_GOlmc_DEGs)
print("length(NKcells_GOlmc_DEGs)")
length(NKcells_GOlmc_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOlmc$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOlmc$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOlmc_breaks <- c(
    (min(NKcells_GOlmc$avg_log2FC) + (min(NKcells_GOlmc$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOlmc$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOlmc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOlmc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOlmc$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOlmc$avg_log2FC) + (max(NKcells_GOlmc$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOlmc_breaks<-round(NKcells_GOlmc_breaks, 2)  # Round normally for the rest
NKcells_GOlmc_breaks<- c(min_value, NKcells_GOlmc_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOlmc_breaks")
NKcells_GOlmc_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOlcm_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOlmc_breaks<-round(NKcells_GOlmc_breaks, 2)  # Round normally for the rest
print("NKcells_GOlmc_breaks")
NKcells_GOlmc_breaks
}
all_Clusters <- unique(NKcells_GOlmc$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOleukocyte_mediated_cytotoxicity_dotPlotCSP.pdf", width=6 ,height=6)
ggplot( NKcells_GOlmc, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOlmc_breaks,
      labels =  NKcells_GOlmc_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte mediated cytotoxicity' in NK cells")  # Add title
dev.off()


print("NKcells lymphocyte proliferation Plot")
lymphocyte_proliferation <- NKcells_all %>% 
  dplyr::filter(Description == "lymphocyte proliferation") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

NKcells_GOlp <- NKcells_DE_all %>%
 dplyr::filter(Human.gene.stable.ID %in% lymphocyte_proliferation) 

write.table(NKcells_GOlp, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOlp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcells_GOlp_DEGs<-NKcells_GOlp$gene
NKcells_GOlp_DEGs <- unique(NKcells_GOlp_DEGs)
print("length(NKcells_GOlp_DEGs)")
length(NKcells_GOlp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOlp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOlp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOlp_breaks <- c(
    (min(NKcells_GOlp$avg_log2FC) + (min(NKcells_GOlp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOlp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOlp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOlp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOlp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOlp$avg_log2FC) + (max(NKcells_GOlp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOlp_breaks<-round(NKcells_GOlp_breaks, 2)  # Round normally for the rest
NKcells_GOlp_breaks<- c(min_value, NKcells_GOlp_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOlp_breaks")
NKcells_GOlp_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOlp_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOlp_breaks<-round(NKcells_GOlp_breaks, 2)  # Round normally for the rest
print("NKcells_GOlp_breaks")
NKcells_GOlp_breaks
}
all_Clusters <- unique(NKcells_GOlp$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOlymphocyte_proliferation_dotPlotCSP.pdf", width=5 ,height=13)
ggplot( NKcells_GOlp, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOlp_breaks,
      labels =  NKcells_GOlp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'lymphocyte proliferation' in NK cells")  # Add title
dev.off()

print("NKcells immune response Plot")
immune_response <- NKcells_all %>% 
  dplyr::filter(Description == "immune response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

NKcells_GOir_all <- NKcells_DE_all %>%
  dplyr::filter(Human.gene.stable.ID %in% immune_response) 

NKcells_GOir <- NKcells_DE_all %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(pct.1 >=60 | pct.2 >= 60, abs(avg_log2FC) >= 2),Human.gene.stable.ID %in% immune_response) %>% 
  dplyr::ungroup()
write.table(NKcells_GOir_all, "/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/NKcells_GOir_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
immune_response_paper <- c("S100A13","RAC2","NKG7","KLRK1","HSPA8","HCST","CORO1A","CD74","CD2","B2M","FCER1G","TRYOBP")
NKcells_GOir_paper <- NKcells_DE_all %>%
  dplyr::filter(gene %in% immune_response_paper)
#make vector of genes
NKcells_GOir_paper_DEGs<-NKcells_GOir_paper$gene
NKcells_GO_DEGs <- unique(NKcells_GOir_paper_DEGs)
print("length(NKcells_GOir_paper_DEGs)")
length(NKcells_GOir_paper_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcells_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcells_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcells_GOir_paper_breaks <- c(
    (min(NKcells_GOir_paper$avg_log2FC) + (min(NKcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcells_GOir_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcells_GOir_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcells_GOir_paper$avg_log2FC) + (max(NKcells_GOir_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcells_GOir_paper_breaks<-round(NKcells_GOir_paper_breaks, 2)  # Round normally for the rest
NKcells_GOir_paper_breaks<- c(min_value, NKcells_GOir_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcells_GOir_paper_breaks")
NKcells_GOir_paper_breaks
}
#if min_value is >= 0, then make breaks this way
if (min_value >= 0) {
midpoint <- (min_value + max_value) / 2 
# Calculate the midpoint between min and the midpoint
midpoint_min <- (min_value + midpoint) / 2
# Calculate the midpoint between max and the midpoint
midpoint_max <- (max_value + midpoint) / 2
# Calculate two additional values between min and midpoint
# Calculate the value between min_value and midpoint_min
min_value2 <- (min_value + midpoint_min) / 2
min_value4 <- (midpoint_min + midpoint) / 2  
# Calculate the value between max_value and midpoint_max
max_value4  <- (midpoint_max + max_value) / 2
# Calculate the value between midpoint_max and midpoint
max_value2 <- (midpoint_max + midpoint) / 2
# Create the breaks vector
NKcells_GOir_paper_breaks <- c(
    min_value,       # Minimum value
    min_value2,     # Value between min and midpoint_min
    midpoint_min,    # Midpoint between min and midpoint
    min_value4,     # Value between midpoint_min and midpoint
    midpoint,        # Midpoint between min and max
    max_value2,     # Value between midpoint and midpoint_max
    midpoint_max,    # Midpoint between max and midpoint
    max_value4,     # Value between midpoint_max and max
    max_value        # Maximum value
)
NKcells_GOir_paper_breaks<-round(NKcells_GOir_paper_breaks, 2)  # Round normally for the rest
print("NKcells_GOir_paper_breaks")
NKcells_GOir_paper_breaks
}
all_Clusters <- unique(NKcells_GOir_paper$Cluster)

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOimmune_response_dotPlotCSP.pdf", width=5 ,height=6.5)
ggplot( NKcells_GOir_paper, aes(x = Cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcells_GOir_paper_breaks,
      labels =  NKcells_GOir_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Clusters ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),  # Increase size and set color to black
    axis.text.y = element_text(face = "italic",size = 12, color = "black"),  # Increase size and set color to black
    axis.title.x = element_text(size = 14, color = "black"),  # Increase size and set color to black
    axis.title.y = element_text(size = 14, color = "black"),  # Increase size and set color to black
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),  # Align title to the middle
    plot.title.position = "plot",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14),# Increase legend title size
    legend.key.size = unit(1, "cm"),  # Increase the size of the color bar  
     panel.spacing = unit(0.01, "lines")  # Reduce space between columns
) +  
    labs(x = "Cluster", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in NK cells")  # Add title
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