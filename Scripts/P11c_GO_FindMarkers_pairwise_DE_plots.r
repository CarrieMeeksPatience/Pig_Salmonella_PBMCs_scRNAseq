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


DE_all_sig <- read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Sal_FindMarkers_pairwise_SIG_DEGs_allCelltypes.txt", sep = "\t", header = TRUE)
DE_all_sig$Comparison <- gsub("D2v0", "2 DPI vs 0 DPI", DE_all_sig$Comparison)
DE_all_sig$Comparison <- gsub("D8v2", "8 DPI vs 2 DPI", DE_all_sig$Comparison)
DE_all_sig$Comparison <- gsub("D8v0", "8 DPI vs 0 DPI", DE_all_sig$Comparison)
DE_all_sig$pct.1 <- DE_all_sig$pct.1 * 100
DE_all_sig$pct.2 <- DE_all_sig$pct.2 * 100
DE_all_sig$"pct.1-pct.2"<- DE_all_sig$pct.1 - DE_all_sig$pct.2
ORG<- read.csv(file ="/Annotation_files/PigToHuman_GeneOrthos_v11_1_97_scGenes.csv", header = T,row.names=1)
colnames(ORG)[11] <- "Gene"
# Replace "_" with "-" in the Gene column
ORG$Gene <- gsub("_", "-", ORG$Gene)
# For duplicates in the Gene column, keep the row with the highest X.id..query.gene.identical.to.target.Human.gene
#ORG2 <- ORG %>% group_by(Gene) %>% filter(X.id..query.gene.identical.to.target.Human.gene == max(X.id..query.gene.identical.to.target.Human.gene)) %>% ungroup()
#subset the columns needed for annotation
ORG2.subset <- ORG %>% dplyr::select(c(1,11,7,8,10))
ORG2.subset$Gene <- as.character(ORG2.subset$Gene)
DE_all_sig<-left_join(DE_all_sig,ORG2.subset,by=c("gene"="Gene"))
DE_all_sig <- unique(DE_all_sig)
# Filter for Monocytes
Monocyte_day <- DE_all_sig %>% filter(Celltype == "Monocytes")
# Set the desired order for the Comparison factor
Monocyte_day$Comparison <- factor(Monocyte_day$Comparison, levels = c("2DPI vs 0DPI", "8DPI vs 2DPI", "8DPI vs 0DPI"))

all_Comparisons <- unique(Monocyte_day$Comparison)



Monocyte_GOresults_df <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", header = TRUE)
Monocyte_GOresults_df$Description <- trimws(as.character(Monocyte_GOresults_df$Description))

print("Monocyte positive regulation of defense response Plot")
positive_regulation_of_defense_response_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "positive regulation of defense response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()

  Monocyte_GOprodr <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% positive_regulation_of_defense_response_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() 

write.table(Monocyte_GOprodr, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOprodr_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOprodr_DEGs<-Monocyte_GOprodr$gene
Monocyte_GOprodr_DEGs <- unique(Monocyte_GOprodr_DEGs)
print("length(Monocyte_GOprodr_DEGs)")
length(Monocyte_GOprodr_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOprodr$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOprodr$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOprodr_breaks <- c(
    (min(Monocyte_GOprodr$avg_log2FC) + (min(Monocyte_GOprodr$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOprodr$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOprodr$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOprodr$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOprodr$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOprodr$avg_log2FC) + (max(Monocyte_GOprodr$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOprodr_breaks<-round(Monocyte_GOprodr_breaks, 2)  # Round normally for the rest
Monocyte_GOprodr_breaks<- c(min_value, Monocyte_GOprodr_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOprodr_breaks")
Monocyte_GOprodr_breaks
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
Monocyte_GOprodr_breaks <- c(
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
Monocyte_GOprodr_breaks<-round(Monocyte_GOprodr_breaks, 2)  # Round normally for the rest
print("Monocyte_GOprodr_breaks")
Monocyte_GOprodr_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOpositive_regulation_of_defense_response_dotPlotCSP.pdf", width=5 ,height=12)
ggplot( Monocyte_GOprodr, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOprodr_breaks,
      labels =  Monocyte_GOprodr_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'positive regulation of defense response' in Monocytes")  # Add title
dev.off()

print("Monocyte regulation of response to stimulus Plot")
regulation_of_response_to_stimulus_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "regulation of response to stimulus") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOrorts <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(pct.1 >= 40 & abs(avg_log2FC) >= 0.5 | pct.2 >= 40 & abs(avg_log2FC) >= 0.5, Human.gene.stable.ID %in% regulation_of_response_to_stimulus_Monocytes) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
table(Monocyte_GOrorts$Comparison)   
write.table(Monocyte_GOrorts, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOrorts_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)      
#make vector of genes
Monocyte_GOrorts_DEGs<-Monocyte_GOrorts$gene
Monocyte_GOrorts_DEGs <- unique(Monocyte_GOrorts_DEGs)
print("length(Monocyte_GOrorts_DEGs)")
length(Monocyte_GOrorts_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOrorts$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOrorts$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOrorts_breaks <- c(
    (min(Monocyte_GOrorts$avg_log2FC) + (min(Monocyte_GOrorts$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOrorts$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOrorts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOrorts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOrorts$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOrorts$avg_log2FC) + (max(Monocyte_GOrorts$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOrorts_breaks<-round(Monocyte_GOrorts_breaks, 2)  # Round normally for the rest
Monocyte_GOrorts_breaks<- c(min_value, Monocyte_GOrorts_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOrorts_breaks")
Monocyte_GOrorts_breaks
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
Monocyte_GOrorts_breaks <- c(
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
Monocyte_GOrorts_breaks<-round(Monocyte_GOrorts_breaks, 2)  # Round normally for the rest
print("Monocyte_GOrorts_breaks")
Monocyte_GOrorts_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOregulation_of_response_to stimulus_dotPlotCSP.pdf", width=5 ,height=13)
ggplot( Monocyte_GOrorts, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOrorts_breaks,
      labels =  Monocyte_GOrorts_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'regulation of response to stimulus' in Monocytes")  # Add title
dev.off()

print("Monocyte immune system process Plot")
immune_system_process_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "immune system process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOisp <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_system_process_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOisp, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOisp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOisp_DEGs<-Monocyte_GOisp$gene
Monocyte_GOisp_DEGs <- unique(Monocyte_GOisp_DEGs)
print("length(Monocyte_GOisp_DEGs)")
length(Monocyte_GOisp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOisp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOisp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOisp_breaks <- c(
    (min(Monocyte_GOisp$avg_log2FC) + (min(Monocyte_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOisp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOisp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOisp$avg_log2FC) + (max(Monocyte_GOisp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOisp_breaks<-round(Monocyte_GOisp_breaks, 2)  # Round normally for the rest
Monocyte_GOisp_breaks<- c(min_value, Monocyte_GOisp_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOisp_breaks")
Monocyte_GOisp_breaks
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
Monocyte_GOisp_breaks <- c(
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
Monocyte_GOisp_breaks<-round(Monocyte_GOisp_breaks, 2)  # Round normally for the rest
print("Monocyte_GOisp_breaks")
Monocyte_GOisp_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOimmune_system_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOisp, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOisp_breaks,
      labels =  Monocyte_GOisp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune system process' in Monocytes")  # Add title
dev.off()

print("Monocyte immune system process paperPlot")
immune_system_process_Monocytes_paper <- c("TAP1","SELL","S100A9","S100A8","S100A12","PSMB8","PSMB10","MEF2C","MAP3K1","LYZ","ITGA4","IFI6","HSPA8","HADC9","GAPDH","CALR","AIF1")
  Monocyte_GOisp_paper <- Monocyte_day %>%
  dplyr::filter(gene %in% immune_system_process_Monocytes_paper)
write.table(Monocyte_GOisp_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOisp_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOisp_paper_DEGs<-Monocyte_GOisp_paper$gene
Monocyte_GOisp_paper_DEGs <- unique(Monocyte_GOisp_paper_DEGs)
print("length(Monocyte_GOisp_paper_DEGs)")
length(Monocyte_GOisp_paper_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOisp_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOisp_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOisp_paper_breaks <- c(
    (min(Monocyte_GOisp_paper$avg_log2FC) + (min(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOisp_paper$avg_log2FC) + (max(Monocyte_GOisp_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOisp_paper_breaks<-round(Monocyte_GOisp_paper_breaks, 2)  # Round normally for the rest
Monocyte_GOisp_paper_breaks<- c(min_value, Monocyte_GOisp_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOisp_paper_breaks")
Monocyte_GOisp_paper_breaks
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
Monocyte_GOisp_paper_breaks <- c(
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
Monocyte_GOisp_paper_breaks<-round(Monocyte_GOisp_paper_breaks, 2)  # Round normally for the rest
print("Monocyte_GOisp_paper_breaks")
Monocyte_GOisp_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOimmune_system_process_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOisp_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOisp_paper_breaks,
      labels =  Monocyte_GOisp_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune system process' in Monocytes for paper")  # Add title
dev.off()

print("Monocyte immune responses Plot")
immune_response_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "immune response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOir <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_response_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOir, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOir_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOir_DEGs<-Monocyte_GOir$gene
Monocyte_GOir_DEGs <- unique(Monocyte_GOir_DEGs)
print("length(Monocyte_GOir_DEGs)")
length(Monocyte_GOir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOir$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOir$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOir_breaks <- c(
    (min(Monocyte_GOir$avg_log2FC) + (min(Monocyte_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOir$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOir$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOir$avg_log2FC) + (max(Monocyte_GOir$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOir_breaks<-round(Monocyte_GOir_breaks, 2)  # Round normally for the rest
Monocyte_GOir_breaks<- c(min_value, Monocyte_GOir_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOir_breaks")
Monocyte_GOir_breaks
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
Monocyte_GOir_breaks <- c(
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
Monocyte_GOir_breaks<-round(Monocyte_GOir_breaks, 2)  # Round normally for the rest
print("Monocyte_GOir_breaks")
Monocyte_GOir_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOimmune_response_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOir, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOir_breaks,
      labels =  Monocyte_GOir_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in Monocytes")  # Add title
dev.off()
print("Monocyte immune responses paper Plot")
immune_response_Monocytes_paper <- c("STAT1","SLA-DQB1","S100A9","S100A8","PTPRC","LYZ","IFITM3","GAPDH","FCN1","CTSS","CTSC","CD74","B2M","MAP3K1")
  Monocyte_GOir_paper<- Monocyte_day %>%
  dplyr::filter(gene %in% immune_response_Monocytes_paper) 
write.table(Monocyte_GOir_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOir_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOir_DEGs<-Monocyte_GOir_paper$gene
Monocyte_GOir_DEGs <- unique(Monocyte_GOir_DEGs)
print("length(Monocyte_GOir_paper_DEGs)")
length(Monocyte_GOir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOir_paper_breaks <- c(
    (min(Monocyte_GOir_paper$avg_log2FC) + (min(Monocyte_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOir_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOir_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOir_paper$avg_log2FC) + (max(Monocyte_GOir_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOir_paper_breaks<-round(Monocyte_GOir_paper_breaks, 2)  # Round normally for the rest
Monocyte_GOir_paper_breaks<- c(min_value, Monocyte_GOir_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOir_paper_breaks")
Monocyte_GOir_paper_breaks
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
Monocyte_GOir_paper_breaks <- c(
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
Monocyte_GOir_paper_breaks<-round(Monocyte_GOir_paper_breaks, 2)  # Round normally for the rest
print("Monocyte_GOir_paper_breaks")
Monocyte_GOir_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOimmune_response_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOir_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOir_paper_breaks,
      labels =  Monocyte_GOir_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in Monocytes for paper")  # Add title
dev.off()


print("Monocyte endocytosis Plot")
endocytosis_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "endocytosis") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOe <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% endocytosis_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOe, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOe_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOe_DEGs<-Monocyte_GOe$gene
Monocyte_GOe_DEGs <- unique(Monocyte_GOe_DEGs)
print("length(Monocyte_GOe_DEGs)")
length(Monocyte_GOe_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOe$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOe$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOe_breaks <- c(
    (min(Monocyte_GOe$avg_log2FC) + (min(Monocyte_GOe$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOe$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOe$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOe$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOe$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOe$avg_log2FC) + (max(Monocyte_GOe$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOe_breaks<-round(Monocyte_GOe_breaks, 2)  # Round normally for the rest
Monocyte_GOe_breaks<- c(min_value, Monocyte_GOe_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOe_breaks")
Monocyte_GOe_breaks
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
Monocyte_GOe_breaks <- c(
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
Monocyte_GOe_breaks<-round(Monocyte_GOe_breaks, 2)  # Round normally for the rest
print("Monocyte_GOe_breaks")
Monocyte_GOe_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOendocytosis_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOe, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOe_breaks,
      labels =  Monocyte_GOe_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'endocytosis' in Monocytes")  # Add title
dev.off()
print("Monocyte endocytosis Plot")
endocytosis_Monocytes_paper <- c("TBC1D5","SH3KBP1","RARA","RACK1","PTPRC","PPP3CA","PIK3CB","LYST","LRRK2","ITGA4","IL10RA","EZR","DPYSL2","DOCK2","CDC42SE2","CD151","BMP2K")
  Monocyte_GOe_paper <- Monocyte_day %>%
  dplyr::filter(gene %in% endocytosis_Monocytes_paper)
write.table(Monocyte_GOe_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOe_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOe_DEGs<-Monocyte_GOe_paper$gene
Monocyte_GOe_DEGs <- unique(Monocyte_GOe_DEGs)
print("length(Monocyte_GOe_paper_DEGs)")
length(Monocyte_GOe_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOe_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOe_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOe_paper_breaks <- c(
    (min(Monocyte_GOe_paper$avg_log2FC) + (min(Monocyte_GOe_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOe_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOe_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOe_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOe_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOe_paper$avg_log2FC) + (max(Monocyte_GOe_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOe_paper_breaks<-round(Monocyte_GOe_paper_breaks, 2)  # Round normally for the rest
Monocyte_GOe_paper_breaks<- c(min_value, Monocyte_GOe_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOe_paper_breaks")
Monocyte_GOe_paper_breaks
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
Monocyte_GOe_paper_breaks <- c(
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
Monocyte_GOe_paper_breaks<-round(Monocyte_GOe_paper_breaks, 2)  # Round normally for the rest
print("Monocyte_GOe_paper_breaks")
Monocyte_GOe_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOendocytosis_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOe_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOe_paper_breaks,
      labels =  Monocyte_GOe_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'endocytosis' in Monocytes for paper")  # Add title
dev.off()


print("Monocyte cell migration Plot")
cell_migration_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "cell migration") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOcm <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% cell_migration_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOcm, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOcm_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOcm_DEGs<-Monocyte_GOcm$gene
Monocyte_GOcm_DEGs <- unique(Monocyte_GOcm_DEGs)
print("length(Monocyte_GOcm_DEGs)")
length(Monocyte_GOcm_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOcm$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOcm$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOcm_breaks <- c(
    (min(Monocyte_GOcm$avg_log2FC) + (min(Monocyte_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOcm$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOcm$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOcm$avg_log2FC) + (max(Monocyte_GOcm$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOcm_breaks<-round(Monocyte_GOcm_breaks, 2)  # Round normally for the rest
Monocyte_GOcm_breaks<- c(min_value, Monocyte_GOcm_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOcm_breaks")
Monocyte_GOcm_breaks
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
Monocyte_GOcm_breaks <- c(
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
Monocyte_GOcm_breaks<-round(Monocyte_GOcm_breaks, 2)  # Round normally for the rest
print("Monocyte_GOcm_breaks")
Monocyte_GOcm_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOcell_migration_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOcm, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOcm_breaks,
      labels =  Monocyte_GOcm_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cell migration' in Monocytes")  # Add title
dev.off()

print("Monocyte cytokine production Plot")
cytokine_production_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "cytokine production") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOcp <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% cytokine_production_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOcp, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOcp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
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

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOcytokine_production_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOcp, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
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

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cytokine production' in Monocytes")  # Add title
dev.off()

print("Monocyte interferon-mediated signaling pathway Plot")
interferon_signaling_Monocytes <- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "interferon-mediated signaling pathway") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOimsp <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% interferon_signaling_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOimsp, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOimsp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOimsp_DEGs<-Monocyte_GOimsp$gene
Monocyte_GOimsp_DEGs <- unique(Monocyte_GOimsp_DEGs)
print("length(Monocyte_GOimsp_DEGs)")
length(Monocyte_GOimsp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOimsp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOimsp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOimsp_breaks <- c(
    (min(Monocyte_GOimsp$avg_log2FC) + (min(Monocyte_GOimsp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOimsp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOimsp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOimsp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOimsp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOimsp$avg_log2FC) + (max(Monocyte_GOimsp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOimsp_breaks<-round(Monocyte_GOimsp_breaks, 2)  # Round normally for the rest
Monocyte_GOimsp_breaks<- c(min_value, Monocyte_GOimsp_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOimsp_breaks")
Monocyte_GOimsp_breaks
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
Monocyte_GOimsp_breaks <- c(
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
Monocyte_GOimsp_breaks<-round(Monocyte_GOimsp_breaks, 2)  # Round normally for the rest
print("Monocyte_GOimsp_breaks")
Monocyte_GOimsp_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOinterferon_mediated_signaling_pathway_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOimsp, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOimsp_breaks,
      labels =  Monocyte_GOimsp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'interferon-mediated signaling pathway' in Monocytes")  # Add title
dev.off()

print("Monocyte positive regulation of cell activation Plot")
positive_regulation_of_cell_activation_Monocytes<- Monocyte_GOresults_df %>% 
  dplyr::filter(Description == "positive regulation of cell activation") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Monocyte_GOproca <- Monocyte_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% positive_regulation_of_cell_activation_Monocytes)) %>%  # Check if any row satisfies the condition
      ungroup() # Ungroup after filtering
write.table(Monocyte_GOproca, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GOproca_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_GOproca_DEGs<-Monocyte_GOproca$gene
Monocyte_GOproca_DEGs <- unique(Monocyte_GOproca_DEGs)
print("length(Monocyte_GOproca_DEGs)")
length(Monocyte_GOproca_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_GOproca$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_GOproca$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_GOproca_breaks <- c(
    (min(Monocyte_GOproca$avg_log2FC) + (min(Monocyte_GOproca$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_GOproca$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_GOproca$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_GOproca$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_GOproca$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_GOproca$avg_log2FC) + (max(Monocyte_GOproca$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_GOproca_breaks<-round(Monocyte_GOproca_breaks, 2)  # Round normally for the rest
Monocyte_GOproca_breaks<- c(min_value, Monocyte_GOproca_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_GOproca_breaks")
Monocyte_GOproca_breaks
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
Monocyte_GOproca_breaks <- c(
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
Monocyte_GOproca_breaks<-round(Monocyte_GOproca_breaks, 2)  # Round normally for the rest
print("Monocyte_GOproca_breaks")
Monocyte_GOproca_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DE_GOpositive_regulation_of_cell_activation_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Monocyte_GOproca, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_GOproca_breaks,
      labels =  Monocyte_GOproca_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'positive regulation of cell activation' in Monocytes")  # Add title
dev.off()

print("Monocyte DE plot")
Monocyte_DEGs <-c("ZEB2","ZDHHC13","STAT1","SLA-DRA","SLA-DB1","S100A9","S100A8","RUNX1","RACK1","PTPRC","LYZ","LGALS1","IFITM3","GAPDH","FLNA","FCN1","CTSS","CTSC","CD74","CD163","CALR","B2M","ACTG1","ACTB")
Monocyte_DE_plot <- Monocyte_day %>%
  dplyr::filter(gene %in% Monocyte_DEGs)
write.table(Monocyte_DE_plot, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Monocyte_GO_DEGs_paper_dotPlot.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Monocyte_DE_plot_DEGs<-Monocyte_DE_plot$gene
Monocyte_DE_plot_DEGs <- unique(Monocyte_DE_plot_DEGs)
print("length(Monocyte_DE_plot_DEGs)")
length(Monocyte_DE_plot_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Monocyte_DE_plot$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Monocyte_DE_plot$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Monocyte_DE_plot_breaks <- c(
    (min(Monocyte_DE_plot$avg_log2FC) + (min(Monocyte_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Monocyte_DE_plot$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Monocyte_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Monocyte_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Monocyte_DE_plot$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Monocyte_DE_plot$avg_log2FC) + (max(Monocyte_DE_plot$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Monocyte_DE_plot_breaks<-round(Monocyte_DE_plot_breaks, 2)  # Round normally for the rest
Monocyte_DE_plot_breaks<- c(min_value, Monocyte_DE_plot_breaks, max_value, 0)  # Add min and max to the breaks    
print("Monocyte_DE_plot_breaks")
Monocyte_DE_plot_breaks
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
Monocyte_DE_plot_breaks <- c(
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
Monocyte_DE_plot_breaks<-round(Monocyte_DE_plot_breaks, 2)  # Round normally for the rest
print("Monocyte_DE_plot_breaks")
Monocyte_DE_plot_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Monocytes_DEGs_paper_dotPlotCSP.pdf", width=5 ,height=13)
ggplot( Monocyte_DE_plot, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Monocyte_DE_plot_breaks,
      labels =  Monocyte_DE_plot_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in Monocytes")  # Add title
dev.off()

# Filter for NK cells
NKcell_day <- DE_all_sig %>% filter(Celltype == "NK cells")
print("table(NKcell_day$Comparison)")
table(NKcell_day$Comparison)
# Set the desired order for the Comparison factor
NKcell_day$Comparison <- factor(NKcell_day$Comparison, levels = c("2DPI vs 0DPI", "8DPI vs 2DPI", "8DPI vs 0DPI"))
# Set the desired order for the Comparison factor
all_Comparisons <- unique(NKcell_day$Comparison)
all_GOresults_df <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_celltypes_pos_enrichGO_all_2026_04_14.txt", sep = "\t", header = TRUE)
# Filter for Monocytes
NKcell_GOresults_df <- DE_all_sig %>% filter(Celltype == "NK cells")

NKcell_GOresults_df <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_NKcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", header = TRUE)
NKcell_GOresults_df$Description <- trimws(as.character(NKcell_GOresults_df$Description))

print("NKcell response to stimulus Plot")
response_to_stimulus_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "response to stimulus") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOrts <- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% response_to_stimulus_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOrts, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOrts_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOrts_DEGs<-NKcell_GOrts$gene
NKcell_GOrts_DEGs <- unique(NKcell_GOrts_DEGs)
print("length(NKcell_GOrts_DEGs)")
length(NKcell_GOrts_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOrts$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOrts$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOrts_breaks <- c(
    (min(NKcell_GOrts$avg_log2FC) + (min(NKcell_GOrts$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOrts$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOrts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOrts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOrts$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOrts$avg_log2FC) + (max(NKcell_GOrts$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOrts_breaks<-round(NKcell_GOrts_breaks, 2)  # Round normally for the rest
NKcell_GOrts_breaks<- c(min_value, NKcell_GOrts_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOrts_breaks")
NKcell_GOrts_breaks
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
NKcell_GOrts_breaks <- c(
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
NKcell_GOrts_breaks<-round(NKcell_GOrts_breaks, 2)  # Round normally for the rest
print("NKcell_GOrts_breaks")
NKcell_GOrts_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOresponse_to_stimulus_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( NKcell_GOrts, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOrts_breaks,
      labels =  NKcell_GOrts_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to stimulus' in NK cells")  # Add title
dev.off()

print("NKcell response to other organism Plot")
response_to_other_organism_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "response to other organism") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOrtoo <- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% response_to_other_organism_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOrtoo, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOrtoo_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOrtoo_DEGs<-NKcell_GOrtoo$gene
NKcell_GOrtoo_DEGs <- unique(NKcell_GOrtoo_DEGs)
print("length(NKcell_GOrtoo_DEGs)")
length(NKcell_GOrtoo_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOrtoo$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOrtoo$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOrtoo_breaks <- c(
    (min(NKcell_GOrtoo$avg_log2FC) + (min(NKcell_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOrtoo$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOrtoo$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOrtoo$avg_log2FC) + (max(NKcell_GOrtoo$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOrtoo_breaks<-round(NKcell_GOrtoo_breaks, 2)  # Round normally for the rest
NKcell_GOrtoo_breaks<- c(min_value, NKcell_GOrtoo_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOrtoo_breaks")
NKcell_GOrtoo_breaks
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
NKcell_GOrtoo_breaks <- c(
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
NKcell_GOrtoo_breaks<-round(NKcell_GOrtoo_breaks, 2)  # Round normally for the rest
print("NKcell_GOrtoo_breaks")
NKcell_GOrtoo_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOresponse_to_other_organism_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( NKcell_GOrtoo, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOrtoo_breaks,
      labels =  NKcell_GOrtoo_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to other organism' in NK cells")  # Add title
dev.off()

print("NKcell response to other organism paper Plot")
response_to_other_organism_NKcells_paper <- c("S100A9","PLAC8","NKG7","LGALS8","IKZF3","GNLY","FCER1G","CORO1A","CFL1","CD74","CCL5","B2M","ACTG1")
  NKcell_GOrtoo_paper <- NKcell_day %>%
  dplyr::filter(gene %in% response_to_other_organism_NKcells_paper)
write.table(NKcell_GOrtoo_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOrtoo_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOrtoo_DEGs<-NKcell_GOrtoo_paper$gene
NKcell_GOrtoo_DEGs <- unique(NKcell_GOrtoo_DEGs)
print("length(NKcell_GOrtoo_DEGs)")
length(NKcell_GOrtoo_DEGs)
#make breaks

# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOrtoo_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOrtoo_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOrtoo_paper_breaks <- c(
    (min(NKcell_GOrtoo_paper$avg_log2FC) + (min(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOrtoo_paper$avg_log2FC) + (max(NKcell_GOrtoo_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOrtoo_paper_breaks<-round(NKcell_GOrtoo_paper_breaks, 2)  # Round normally for the rest
NKcell_GOrtoo_paper_breaks<- c(min_value, NKcell_GOrtoo_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOrtoo_paper_breaks")
NKcell_GOrtoo_paper_breaks
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
NKcell_GOrtoo_paper_breaks <- c(
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
NKcell_GOrtoo_paper_breaks<-round(NKcell_GOrtoo_paper_breaks, 2)  # Round normally for the rest
print("NKcell_GOrtoo_paper_breaks")
NKcell_GOrtoo_paper_breaks
}
all_Comparisons_NKcell_GOrtoo_paper <- unique(NKcell_GOrtoo_paper$Comparison)
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOresponse_to_other_organism_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( NKcell_GOrtoo_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOrtoo_paper_breaks,
      labels =  NKcell_GOrtoo_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons_NKcell_GOrtoo_paper ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to other organism' in NK cells for paper")  # Add title
dev.off()

print("NKcell response to biotic stimulus Plot")
response_to_biotic_stimulus_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "response to biotic stimulus") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOrtbs <- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% response_to_biotic_stimulus_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOrtbs, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOrtbs_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOrtbs_DEGs<-NKcell_GOrtbs$gene
NKcell_GOrtbs_DEGs <- unique(NKcell_GOrtbs_DEGs)
print("length(NKcell_GOrtbs_DEGs)")
length(NKcell_GOrtbs_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOrtbs$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOrtbs$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOrtbs_breaks <- c(
    (min(NKcell_GOrtbs$avg_log2FC) + (min(NKcell_GOrtbs$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOrtbs$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOrtbs$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOrtbs$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOrtbs$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOrtbs$avg_log2FC) + (max(NKcell_GOrtbs$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOrtbs_breaks<-round(NKcell_GOrtbs_breaks, 2)  # Round normally for the rest
NKcell_GOrtbs_breaks<- c(min_value, NKcell_GOrtbs_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOrtbs_breaks")
NKcell_GOrtbs_breaks
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
NKcell_GOrtbs_breaks <- c(
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
NKcell_GOrtbs_breaks<-round(NKcell_GOrtbs_breaks, 2)  # Round normally for the rest
print("NKcell_GOrtbs_breaks")
NKcell_GOrtbs_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOresponse_to_biotic_stimulus_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( NKcell_GOrtbs, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOrtbs_breaks,
      labels =  NKcell_GOrtbs_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to biotic stimulus' in NK cells")  # Add title
dev.off()

print("NKcell leukocyte activation Plot")
leukocyte_activation_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "leukocyte activation") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOla<- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% leukocyte_activation_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOla, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOla_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOla_DEGs<-NKcell_GOla$gene
NKcell_GOla_DEGs <- unique(NKcell_GOla_DEGs)
print("length(NKcell_GOla_DEGs)")
length(NKcell_GOla_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOla$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOla$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOla_breaks <- c(
    (min(NKcell_GOla$avg_log2FC) + (min(NKcell_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOla$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOla$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOla$avg_log2FC) + (max(NKcell_GOla$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOla_breaks<-round(NKcell_GOla_breaks, 2)  # Round normally for the rest
NKcell_GOla_breaks<- c(min_value, NKcell_GOla_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOla_breaks")
NKcell_GOla_breaks
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
NKcell_GOla_breaks <- c(
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
NKcell_GOla_breaks<-round(NKcell_GOla_breaks, 2)  # Round normally for the rest
print("NKcell_GOla_breaks")
NKcell_GOla_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOleukocyte_activation_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(NKcell_GOla, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOla_breaks,
      labels =  NKcell_GOla_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte activation' in NK cells")  # Add title
dev.off()

print("NKcell leukocyte activation paper Plot")
leukocyte_activation_NKcells_paper <- c("TOX","SYK","RUNX1","RORA","PTPRC","PRKCB","PPP3CA","PLCL2","ITGA4","IKZF3","HDAC9","DOCK10")
  NKcell_GOla_paper<- NKcell_day %>%
  dplyr::filter(gene %in% leukocyte_activation_NKcells_paper)
write.table(NKcell_GOla_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOla_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOla_DEGs<-NKcell_GOla_paper$gene
NKcell_GOla_DEGs <- unique(NKcell_GOla_DEGs)
print("length(NKcell_GOla_DEGs)")
length(NKcell_GOla_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOla_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOla_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
  NKcell_GOla_paper_breaks <- c(
    (min(NKcell_GOla_paper$avg_log2FC) + (min(NKcell_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative midpoint and min
    (min(NKcell_GOla_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOla_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOla_paper$avg_log2FC) + (max(NKcell_GOla_paper$avg_log2FC) + 0) / 2) / 2  # Break between positive midpoint and max
  )
  NKcell_GOla_paper_breaks <- round(NKcell_GOla_paper_breaks, 2)  # Round normally for the rest
  NKcell_GOla_paper_breaks <- c(min_value, NKcell_GOla_paper_breaks, max_value, 0)  # Add min and max to the breaks    
  print("NKcell_GOla_paper_breaks")
  NKcell_GOla_paper_breaks
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
NKcell_GOla_paper_breaks <- c(
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
NKcell_GOla_paper_breaks<-round(NKcell_GOla_paper_breaks, 2)  # Round normally for the rest
print("NKcell_GOla_paper_breaks")
NKcell_GOla_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOleukocyte_activation_paper_dotPlotCSP.pdf", width=5 ,height=9)
ggplot(NKcell_GOla_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOla_paper_breaks,
      labels =  NKcell_GOla_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte activation' in NK cells for paper")  # Add title
dev.off()

print("NKcell immune system process Plot")
immune_system_process_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "immune system process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOisp<- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_system_process_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOisp, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOisp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOisp_DEGs<-NKcell_GOisp$gene
NKcell_GOisp_DEGs <- unique(NKcell_GOisp_DEGs)
print("length(NKcell_GOisp_DEGs)")
length(NKcell_GOisp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOisp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOisp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOisp_breaks <- c(
    (min(NKcell_GOisp$avg_log2FC) + (min(NKcell_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOisp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOisp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOisp$avg_log2FC) + (max(NKcell_GOisp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOisp_breaks<-round(NKcell_GOisp_breaks, 2)  # Round normally for the rest
NKcell_GOisp_breaks<- c(min_value, NKcell_GOisp_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOisp_breaks")
NKcell_GOisp_breaks
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
NKcell_GOisp_breaks <- c(
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
NKcell_GOisp_breaks<-round(NKcell_GOisp_breaks, 2)  # Round normally for the rest
print("NKcell_GOisp_breaks")
NKcell_GOisp_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOimmune_system_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(NKcell_GOisp, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOisp_breaks,
      labels =  NKcell_GOisp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune system process' in NK cells")  # Add title
dev.off()

print("NKcell immune response Plot")
immune_response_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "immune response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOir<- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_response_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOir, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOir_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOir_DEGs<-NKcell_GOir$gene
NKcell_GOir_DEGs <- unique(NKcell_GOir_DEGs)
print("length(NKcell_GOir_DEGs)")
length(NKcell_GOir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOir$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOir$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOir_breaks <- c(
    (min(NKcell_GOir$avg_log2FC) + (min(NKcell_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOir$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOir$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOir$avg_log2FC) + (max(NKcell_GOir$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOir_breaks<-round(NKcell_GOir_breaks, 2)  # Round normally for the rest
NKcell_GOir_breaks<- c(min_value, NKcell_GOir_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOir_breaks")
NKcell_GOir_breaks
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
NKcell_GOir_breaks <- c(
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
NKcell_GOir_breaks<-round(NKcell_GOir_breaks, 2)  # Round normally for the rest
print("NKcell_GOir_breaks")
NKcell_GOir_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOimmune_response_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(NKcell_GOir, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOir_breaks,
      labels =  NKcell_GOir_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in NK cells")  # Add title
dev.off()

print("NKcell immune response paper Plot")
immune_response_NKcells_paper <- c("RAP1GAP2","RAC2","PRKCB","ITGA4","GNLY","FCER1G","ETS1","DOCK10","CTSW","CORO1A","CCL5","B2M","ACTG1")
  NKcell_GOir_paper<- NKcell_day %>%
  dplyr::filter(gene %in% immune_response_NKcells_paper)
write.table(NKcell_GOir_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOir_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOir_DEGs<-NKcell_GOir_paper$gene
NKcell_GOir_DEGs <- unique(NKcell_GOir_DEGs)
print("length(NKcell_GOir_DEGs)")
length(NKcell_GOir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOir_paper_breaks <- c(
    (min(NKcell_GOir_paper$avg_log2FC) + (min(NKcell_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOir_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOir_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOir_paper$avg_log2FC) + (max(NKcell_GOir_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOir_paper_breaks<-round(NKcell_GOir_paper_breaks, 2)  # Round normally for the rest
NKcell_GOir_paper_breaks<- c(min_value, NKcell_GOir_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOir_paper_breaks")
NKcell_GOir_paper_breaks
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
NKcell_GOir_paper_breaks <- c(
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
NKcell_GOir_paper_breaks<-round(NKcell_GOir_paper_breaks, 2)  # Round normally for the rest
print("NKcell_GOir_paper_breaks")
NKcell_GOir_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOimmune_response_paper_dotPlotCSP.pdf", width=5 ,height=9)
ggplot(NKcell_GOir_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOir_paper_breaks,
      labels =  NKcell_GOir_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in NK cells for paper")  # Add title
dev.off()

print("NKcell cell migration Plot")
cell_migration_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "cell migration") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOcm_all<- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% cell_migration_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOcm_all, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOcm_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
NKcell_GOcm <- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(pct.1 >= 30 | pct.2 >= 30, abs(avg_log2FC) >= 0.45),Human.gene.stable.ID %in% cell_migration_NKcells) %>%  # Check if any row satisfies the condition
      ungroup() 

#make vector of genes
NKcell_GOcm_DEGs<-NKcell_GOcm$gene
NKcell_GOcm_DEGs <- unique(NKcell_GOcm_DEGs)
print("length(NKcell_GOcm_DEGs)")
length(NKcell_GOcm_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOcm$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOcm$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOcm_breaks <- c(
    (min(NKcell_GOcm$avg_log2FC) + (min(NKcell_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOcm$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOcm$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOcm$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOcm$avg_log2FC) + (max(NKcell_GOcm$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOcm_breaks<-round(NKcell_GOcm_breaks, 2)  # Round normally for the rest
NKcell_GOcm_breaks<- c(min_value, NKcell_GOcm_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOcm_breaks")
NKcell_GOcm_breaks
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
NKcell_GOcm_breaks <- c(
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
NKcell_GOcm_breaks<-round(NKcell_GOcm_breaks, 2)  # Round normally for the rest
print("NKcell_GOcm_breaks")
NKcell_GOcm_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOcell_migration_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(NKcell_GOcm, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOcm_breaks,
      labels =  NKcell_GOcm_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cell migration' in NK cells")  # Add title
dev.off()

print("NKcell apoptotic process Plot")
apoptotic_process_NKcells <- NKcell_GOresults_df %>% 
  dplyr::filter(Description == "apoptotic process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  NKcell_GOap<- NKcell_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% apoptotic_process_NKcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(NKcell_GOap, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GOap_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_GOap_DEGs<-NKcell_GOap$gene
NKcell_GOap_DEGs <- unique(NKcell_GOap_DEGs)
print("length(NKcell_GOap_DEGs)")
length(NKcell_GOap_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_GOap$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_GOap$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_GOap_breaks <- c(
    (min(NKcell_GOap$avg_log2FC) + (min(NKcell_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_GOap$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_GOap$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_GOap$avg_log2FC) + (max(NKcell_GOap$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_GOap_breaks<-round(NKcell_GOap_breaks, 2)  # Round normally for the rest
NKcell_GOap_breaks<- c(min_value, NKcell_GOap_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_GOap_breaks")
NKcell_GOap_breaks
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
NKcell_GOap_breaks <- c(
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
NKcell_GOap_breaks<-round(NKcell_GOap_breaks, 2)  # Round normally for the rest
print("NKcell_GOap_breaks")
NKcell_GOap_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DE_GOapoptotic_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(NKcell_GOap, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_GOap_breaks,
      labels =  NKcell_GOap_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'apoptotic process' in NK cells")  # Add title
dev.off() 

print("NK cells DE plot")
NKcell_DEGs <-c("RUNX1","RAP1GAP2","RACK1","PTPRC","PRKCB","OSBPL8","LGALS8","IKZF3","GNLY","DOCK10","CTSW","CCL5","B2M","ACTG1")
NKcell_DE_plot <- NKcell_day %>%
  dplyr::filter(gene %in% NKcell_DEGs)
write.table(NKcell_DE_plot, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/NKcell_GO_DEGs_paper_dotPlot.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
NKcell_DE_plot_DEGs<-NKcell_DE_plot$gene
NKcell_DE_plot_DEGs <- unique(NKcell_DE_plot_DEGs)
print("length(NKcell_DE_plot_DEGs)")
length(NKcell_DE_plot_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(NKcell_DE_plot$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(NKcell_DE_plot$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
NKcell_DE_plot_breaks <- c(
    (min(NKcell_DE_plot$avg_log2FC) + (min(NKcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(NKcell_DE_plot$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(NKcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(NKcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(NKcell_DE_plot$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(NKcell_DE_plot$avg_log2FC) + (max(NKcell_DE_plot$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
NKcell_DE_plot_breaks<-round(NKcell_DE_plot_breaks, 2)  # Round normally for the rest
NKcell_DE_plot_breaks<- c(min_value, NKcell_DE_plot_breaks, max_value, 0)  # Add min and max to the breaks    
print("NKcell_DE_plot_breaks")
NKcell_DE_plot_breaks
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
NKcell_DE_plot_breaks <- c(
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
NKcell_DE_plot_breaks<-round(NKcell_DE_plot_breaks, 2)  # Round normally for the rest
print("NKcell_DE_plot_breaks")
NKcell_DE_plot_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcells_DEGs_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( NKcell_DE_plot, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  NKcell_DE_plot_breaks,
      labels =  NKcell_DE_plot_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in NK cells")  # Add title
dev.off()

# Filter for B cells
Bcells_day <- DE_all_sig %>% filter(Celltype == "B cells")
# Set the desired order for the Comparison factor
Bcells_day$Comparison <- factor(Bcells_day$Comparison, levels = c("2DPI vs 0DPI", "8DPI vs 2DPI", "8DPI vs 0DPI"))
# Set the desired order for the Comparison factor
all_Comparisons <- unique(Bcells_day$Comparison)

Bcells_GOresults_df <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Bcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", header = TRUE)
Bcells_GOresults_df$Description <- trimws(as.character(Bcells_GOresults_df$Description))


print("Bcells response to other organism Plot")
response_to_other_organism_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "response to other organism") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOrtoo <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% response_to_other_organism_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOrtoo, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOrtoo_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOrtoo_DEGs<-Bcells_GOrtoo$gene
Bcells_GOrtoo_DEGs <- unique(Bcells_GOrtoo_DEGs)
print("length(Bcells_GOrtoo_DEGs)")
length(Bcells_GOrtoo_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOrtoo$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOrtoo$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOrtoo_breaks <- c(
    (min(Bcells_GOrtoo$avg_log2FC) + (min(Bcells_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOrtoo$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOrtoo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOrtoo$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOrtoo$avg_log2FC) + (max(Bcells_GOrtoo$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOrtoo_breaks<-round(Bcells_GOrtoo_breaks, 2)  # Round normally for the rest
Bcells_GOrtoo_breaks<- c(min_value, Bcells_GOrtoo_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOrtoo_breaks")
Bcells_GOrtoo_breaks
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
Bcells_GOrtoo_breaks <- c(
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
Bcells_GOrtoo_breaks<-round(Bcells_GOrtoo_breaks, 2)  # Round normally for the rest
print("Bcells_GOrtoo_breaks")
Bcells_GOrtoo_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOresponse_to_other_organism_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Bcells_GOrtoo, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOrtoo_breaks,
      labels =  Bcells_GOrtoo_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'response to other organism' in B cells")  # Add title
dev.off()

print("Bcells positive regulation of response to stimulus Plot")
positive_regulation_of_response_to_stimulus_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "positive regulation of response to stimulus") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOprorts <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% positive_regulation_of_response_to_stimulus_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOprorts, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOprorts_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOprorts_DEGs<-Bcells_GOprorts$gene
Bcells_GOprorts_DEGs <- unique(Bcells_GOprorts_DEGs)
print("length(Bcells_GOprorts_DEGs)")
length(Bcells_GOprorts_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOprorts$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOprorts$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOprorts_breaks <- c(
    (min(Bcells_GOprorts$avg_log2FC) + (min(Bcells_GOprorts$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOprorts$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOprorts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOprorts$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOprorts$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOprorts$avg_log2FC) + (max(Bcells_GOprorts$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOprorts_breaks<-round(Bcells_GOprorts_breaks, 2)  # Round normally for the rest
Bcells_GOprorts_breaks<- c(min_value, Bcells_GOprorts_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOprorts_breaks")
Bcells_GOprorts_breaks
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
Bcells_GOprorts_breaks <- c(
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
Bcells_GOprorts_breaks<-round(Bcells_GOprorts_breaks, 2)  # Round normally for the rest
print("Bcells_GOprorts_breaks")
Bcells_GOprorts_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOpositive_regulation_of_response_to_stimulus_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Bcells_GOprorts, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOprorts_breaks,
      labels =  Bcells_GOprorts_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'positive regulation of response to stimulus' in B cells")  # Add title
dev.off()

print("Bcells positive regulation of immune response Plot")
positive_regulation_of_immune_response_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "positive regulation of immune response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOproir <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% positive_regulation_of_immune_response_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOproir, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOproir_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOproir_DEGs<-Bcells_GOproir$gene
Bcells_GOproir_DEGs <- unique(Bcells_GOproir_DEGs)
print("length(Bcells_GOproir_DEGs)")
length(Bcells_GOproir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOproir$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOproir$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOproir_breaks <- c(
    (min(Bcells_GOproir$avg_log2FC) + (min(Bcells_GOproir$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOproir$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOproir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOproir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOproir$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOproir$avg_log2FC) + (max(Bcells_GOproir$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOproir_breaks<-round(Bcells_GOproir_breaks, 2)  # Round normally for the rest
Bcells_GOproir_breaks<- c(min_value, Bcells_GOproir_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOproir_breaks")
Bcells_GOproir_breaks
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
Bcells_GOproir_breaks <- c(
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
Bcells_GOproir_breaks<-round(Bcells_GOproir_breaks, 2)  # Round normally for the rest
print("Bcells_GOproir_breaks")
Bcells_GOproir_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOpositive_regulation_of_immune_response_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Bcells_GOproir, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOproir_breaks,
      labels =  Bcells_GOproir_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'positive regulation of immune response' in B cells")  # Add title
dev.off()

print("Bcells leukocyte activation Plot")
leukocyte_activation_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "leukocyte activation") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOla <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% leukocyte_activation_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOla, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOla_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOla_DEGs<-Bcells_GOla$gene
Bcells_GOla_DEGs <- unique(Bcells_GOla_DEGs)
print("length(Bcells_GOla_DEGs)")
length(Bcells_GOla_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOla$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOla$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOla_breaks <- c(
    (min(Bcells_GOla$avg_log2FC) + (min(Bcells_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOla$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOla$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOla$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOla$avg_log2FC) + (max(Bcells_GOla$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOla_breaks<-round(Bcells_GOla_breaks, 2)  # Round normally for the rest
Bcells_GOla_breaks<- c(min_value, Bcells_GOla_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOla_breaks")
Bcells_GOla_breaks
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
Bcells_GOla_breaks <- c(
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
Bcells_GOla_breaks<-round(Bcells_GOla_breaks, 2)  # Round normally for the rest
print("Bcells_GOla_breaks")
Bcells_GOla_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOleukocyte_activation_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( Bcells_GOla, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOla_breaks,
      labels =  Bcells_GOla_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte activation' in B cells")  # Add title
dev.off()

print("Bcells leukocyte activation paper Plot")
leukocyte_activation_Bcells_paper <- c("SLA-DRA","SLA-DQB1","MS4A1","MIF","IKZF3","HDAC9","FAM49B","CD79B","CD79A","CD74","CD151","B2M","ACTB")
  Bcells_GOla_paper <- Bcells_day %>%
  dplyr::filter(gene %in% leukocyte_activation_Bcells_paper)
write.table(Bcells_GOla_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOla_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOla_paper_DEGs<-Bcells_GOla_paper$gene
Bcells_GOla_paper_DEGs <- unique(Bcells_GOla_paper_DEGs)
print("length(Bcells_GOla_paper_DEGs)")
length(Bcells_GOla_paper_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOla_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOla_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOla_paper_breaks <- c(
    (min(Bcells_GOla_paper$avg_log2FC) + (min(Bcells_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOla_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOla_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOla_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOla_paper$avg_log2FC) + (max(Bcells_GOla_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOla_paper_breaks<-round(Bcells_GOla_paper_breaks, 2)  # Round normally for the rest
Bcells_GOla_paper_breaks<- c(min_value, Bcells_GOla_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOla_paper_breaks")
Bcells_GOla_paper_breaks
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
Bcells_GOla_paper_breaks <- c(
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
Bcells_GOla_paper_breaks<-round(Bcells_GOla_paper_breaks, 2)  # Round normally for the rest
print("Bcells_GOla_paper_breaks")
Bcells_GOla_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOleukocyte_activation_paper_dotPlotCSP.pdf", width=5 ,height=9)
ggplot( Bcells_GOla_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOla_paper_breaks,
      labels =  Bcells_GOla_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte activation' in B cells for paper")  # Add title
dev.off()


print("Bcells intracellular signal transduction Plot")
intracellular_signal_transduction_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "intracellular signal transduction") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOist <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% intracellular_signal_transduction_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOist, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOist_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOist_DEGs<-Bcells_GOist$gene
Bcells_GOist_DEGs <- unique(Bcells_GOist_DEGs)
print("length(Bcells_GOist_DEGs)")
length(Bcells_GOist_DEGs)
print("Bcells intracellular signal transduction")
Bcells_GOist_DEGs
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOist$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOist$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOist_breaks <- c(
    (min(Bcells_GOist$avg_log2FC) + (min(Bcells_GOist$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOist$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOist$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOist$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOist$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOist$avg_log2FC) + (max(Bcells_GOist$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOist_breaks<-round(Bcells_GOist_breaks, 2)  # Round normally for the rest
Bcells_GOist_breaks<- c(min_value, Bcells_GOist_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOist_breaks")
Bcells_GOist_breaks
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
Bcells_GOist_breaks <- c(
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
Bcells_GOist_breaks<-round(Bcells_GOist_breaks, 2)  # Round normally for the rest
print("Bcells_GOist_breaks")
Bcells_GOist_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOintracellular_signal_transduction_dotPlotCSP.pdf", width=6 ,height=13)
ggplot(Bcells_GOist, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOist_breaks,
      labels =  Bcells_GOist_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'intracellular signal transduction' in B cells")  # Add title
dev.off()

print("Bcells intracellular signal transduction paper Plot")
intracellular_signal_transduction_Bcells_paper <- c("EZR","VIM","FOXO1","LNPEP","PAG1","LIMD1","MAP3K1","GPR174","ARF6","MCL1","TFRC","RAP1B","PIK3R1","PIK3CA","ZFP36L1","IQGAP1","KRAS")
  Bcells_GOist_paper <- Bcells_day %>%
  dplyr::filter(gene %in% intracellular_signal_transduction_Bcells_paper)
#make vector of genes
Bcells_GOist_DEGs<-Bcells_GOist_paper$gene
Bcells_GOist_DEGs <- unique(Bcells_GOist_DEGs)
print("length(Bcells_GOist_DEGs)")
length(Bcells_GOist_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOist_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOist_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOist_paper_breaks <- c(
    (min(Bcells_GOist_paper$avg_log2FC) + (min(Bcells_GOist_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOist_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOist_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOist_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOist_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOist_paper$avg_log2FC) + (max(Bcells_GOist_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOist_paper_breaks<-round(Bcells_GOist_paper_breaks, 2)  # Round normally for the rest
Bcells_GOist_paper_breaks<- c(min_value, Bcells_GOist_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOist_paper_breaks")
Bcells_GOist_paper_breaks
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
Bcells_GOist_paper_breaks <- c(
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
Bcells_GOist_paper_breaks<-round(Bcells_GOist_paper_breaks, 2)  # Round normally for the rest
print("Bcells_GOist_paper_breaks")
Bcells_GOist_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOintracellular_signal_transduction_paper_dotPlotCSP.pdf", width=6 ,height=13)
ggplot(Bcells_GOist_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOist_paper_breaks,
      labels =  Bcells_GOist_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'intracellular signal transduction' in B cells for paper")  # Add title
dev.off()

print("Bcells immune system process Plot")
immune_system_process_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "immune system process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOisp <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_system_process_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOisp, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOisp_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOisp_DEGs<-Bcells_GOisp$gene
Bcells_GOisp_DEGs <- unique(Bcells_GOisp_DEGs)
print("length(Bcells_GOisp_DEGs)")
length(Bcells_GOisp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOisp$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOisp$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOisp_breaks <- c(
    (min(Bcells_GOisp$avg_log2FC) + (min(Bcells_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOisp$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOisp$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOisp$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOisp$avg_log2FC) + (max(Bcells_GOisp$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOisp_breaks<-round(Bcells_GOisp_breaks, 2)  # Round normally for the rest
Bcells_GOisp_breaks<- c(min_value, Bcells_GOisp_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOisp_breaks")
Bcells_GOisp_breaks
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
Bcells_GOisp_breaks <- c(
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
Bcells_GOisp_breaks<-round(Bcells_GOisp_breaks, 2)  # Round normally for the rest
print("Bcells_GOisp_breaks")
Bcells_GOisp_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOimmune_system_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(Bcells_GOisp, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOisp_breaks,
      labels =  Bcells_GOisp_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune system process' in B cells")  # Add title
dev.off()

print("Bcells immune system process paper Plot")
immune_system_process_Bcells_paper <- c("SLA-DRA","SLA-DQB1","PCBP2","MS4A1","MIF","LY86","KYNU","KLF2","IKZF2","IKZF3","HDAC9","FAM49B","EZR","CIITA","CD79B","CD79A","CD74","CD151","B2M","ACTB")
  Bcells_GOisp_paper <- Bcells_day %>%
  dplyr::filter(gene %in% immune_system_process_Bcells_paper)
write.table(Bcells_GOisp_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOisp_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOisp_DEGs<-Bcells_GOisp_paper$gene
Bcells_GOisp_DEGs <- unique(Bcells_GOisp_DEGs)
print("length(Bcells_GOisp_paper_DEGs)")
length(Bcells_GOisp_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOisp_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOisp_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOisp_paper_breaks <- c(
    (min(Bcells_GOisp_paper$avg_log2FC) + (min(Bcells_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOisp_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOisp_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOisp_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOisp_paper$avg_log2FC) + (max(Bcells_GOisp_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOisp_paper_breaks<-round(Bcells_GOisp_paper_breaks, 2)  # Round normally for the rest
Bcells_GOisp_paper_breaks<- c(min_value, Bcells_GOisp_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOisp_paper_breaks")
Bcells_GOisp_paper_breaks
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
Bcells_GOisp_paper_breaks <- c(
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
Bcells_GOisp_paper_breaks<-round(Bcells_GOisp_paper_breaks, 2)  # Round normally for the rest
print("Bcells_GOisp_paper_breaks")
Bcells_GOisp_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOimmune_system_process_paper_dotPlotCSP.pdf", width=5 ,height=12)
ggplot(Bcells_GOisp_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOisp_paper_breaks,
      labels =  Bcells_GOisp_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune system process' in B cells for paper")  # Add title
dev.off()

print("Bcells immune response Plot")
immune_response_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "immune response") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOir <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% immune_response_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOir, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOir_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOir_DEGs<-Bcells_GOir$gene
Bcells_GOir_DEGs <- unique(Bcells_GOir_DEGs)
print("length(Bcells_GOir_DEGs)")
length(Bcells_GOir_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOir$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOir$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOir_breaks <- c(
    (min(Bcells_GOir$avg_log2FC) + (min(Bcells_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOir$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOir$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOir$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOir$avg_log2FC) + (max(Bcells_GOir$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOir_breaks<-round(Bcells_GOir_breaks, 2)  # Round normally for the rest
Bcells_GOir_breaks<- c(min_value, Bcells_GOir_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOir_breaks")
Bcells_GOir_breaks
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
Bcells_GOir_breaks <- c(
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
Bcells_GOir_breaks<-round(Bcells_GOir_breaks, 2)  # Round normally for the rest
print("Bcells_GOir_breaks")
Bcells_GOir_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOimmune_response_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(Bcells_GOir, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOir_breaks,
      labels =  Bcells_GOir_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in B cells")  # Add title
dev.off()

print("Bcells immune response paperPlot")
immune_response_Bcells_paper <- c("SLA-DRA","SLA-DQB1","PCBP2","MS4A1","MIF","LY86","KYNU","FAM49B","EZR","CIITA","CD79B","CD79A","CD74","B2M")
  Bcells_GOir_paper <- Bcells_day %>%
  dplyr::filter(gene %in% immune_response_Bcells_paper)
write.table(Bcells_GOir_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOir_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOir_paper_DEGs<-Bcells_GOir_paper$gene
Bcells_GOir_paper_DEGs <- unique(Bcells_GOir_paper_DEGs)
print("length(Bcells_GOir_paper_DEGs)")
length(Bcells_GOir_paper_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOir_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOir_paper_breaks <- c(
    (min(Bcells_GOir_paper$avg_log2FC) + (min(Bcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOir_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOir_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOir_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOir_paper$avg_log2FC) + (max(Bcells_GOir_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOir_paper_breaks<-round(Bcells_GOir_paper_breaks, 2)  # Round normally for the rest
Bcells_GOir_paper_breaks<- c(min_value, Bcells_GOir_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOir_paper_breaks")
Bcells_GOir_paper_breaks
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
Bcells_GOir_paper_breaks <- c(
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
Bcells_GOir_paper_breaks<-round(Bcells_GOir_paper_breaks, 2)  # Round normally for the rest
print("Bcells_GOir_paper_breaks")
Bcells_GOir_paper_breaks
}
all_Comparisons<-unique(Bcells_GOir_paper$Comparison)
pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOimmune_response_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(Bcells_GOir_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOir_paper_breaks,
      labels =  Bcells_GOir_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'immune response' in B cells for paper")  # Add title
dev.off()

print("Bcells apoptotic process Plot")
apoptotic_process_Bcells <- Bcells_GOresults_df %>% 
  dplyr::filter(Description == "apoptotic process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  Bcells_GOap <- Bcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% apoptotic_process_Bcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(Bcells_GOap, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcells_GOap_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
Bcells_GOap_DEGs<-Bcells_GOap$gene
Bcells_GOap_DEGs <- unique(Bcells_GOap_DEGs)
print("length(Bcells_GOap_DEGs)")
length(Bcells_GOap_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcells_GOap$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcells_GOap$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcells_GOap_breaks <- c(
    (min(Bcells_GOap$avg_log2FC) + (min(Bcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcells_GOap$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcells_GOap$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcells_GOap$avg_log2FC) + (max(Bcells_GOap$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcells_GOap_breaks<-round(Bcells_GOap_breaks, 2)  # Round normally for the rest
Bcells_GOap_breaks<- c(min_value, Bcells_GOap_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcells_GOap_breaks")
Bcells_GOap_breaks
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
Bcells_GOap_breaks <- c(
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
Bcells_GOap_breaks<-round(Bcells_GOap_breaks, 2)  # Round normally for the rest
print("Bcells_GOap_breaks")
Bcells_GOap_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DE_GOapoptotic_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot(Bcells_GOap, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcells_GOap_breaks,
      labels =  Bcells_GOap_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'apoptotic process' in B cells")  # Add title
dev.off()

print("B cells DE plot")
Bcell_DEGs <-c("VIM","THY1","TFRC","STAT2","STAT1","S100A9","RBPJ","PTPRJ","PSMB10","JAK2","ITGB2","ITGAM","FAM49B","CNR2","CIITA","CD79B","CD79A","CD74","BIRC3","B2M")
Bcell_DE_plot <- Bcells_day %>%
  dplyr::filter(gene %in% Bcell_DEGs)
write.table(Bcell_DE_plot, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/Bcell_GO_DEGs_paper_dotPlot.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#make vector of genes
Bcell_DE_plot_DEGs<-Bcell_DE_plot$gene
Bcell_DE_plot_DEGs <- unique(Bcell_DE_plot_DEGs)
print("length(Bcell_DE_plot_DEGs)")
length(Bcell_DE_plot_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(Bcell_DE_plot$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(Bcell_DE_plot$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
Bcell_DE_plot_breaks <- c(
    (min(Bcell_DE_plot$avg_log2FC) + (min(Bcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(Bcell_DE_plot$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(Bcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(Bcell_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(Bcell_DE_plot$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(Bcell_DE_plot$avg_log2FC) + (max(Bcell_DE_plot$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
Bcell_DE_plot_breaks<-round(Bcell_DE_plot_breaks, 2)  # Round normally for the rest
Bcell_DE_plot_breaks<- c(min_value, Bcell_DE_plot_breaks, max_value, 0)  # Add min and max to the breaks    
print("Bcell_DE_plot_breaks")
Bcell_DE_plot_breaks
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
Bcell_DE_plot_breaks <- c(
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
Bcell_DE_plot_breaks<-round(Bcell_DE_plot_breaks, 2)  # Round normally for the rest
print("Bcell_DE_plot_breaks")
Bcell_DE_plot_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Bcells_DEGs_paper_dotPlotCSP.pdf", width=5 ,height=9)
ggplot( Bcell_DE_plot, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  Bcell_DE_plot_breaks,
      labels =  Bcell_DE_plot_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in B cells")  # Add title
dev.off()

# Filter for CD2- gd T cells
gdTcells_day <- DE_all_sig %>% filter(Celltype == "CD2- gd T cells")
print("dim(gdTcells_day)")
dim(gdTcells_day)
# Set the desired order for the Comparison factor
gdTcells_day$Comparison <- factor(gdTcells_day$Comparison, levels = c("2DPI vs 0DPI", "8DPI vs 2DPI", "8DPI vs 0DPI"))
# Set the desired order for the Comparison factor
all_Comparisons <- unique(gdTcells_day$Comparison)


gdTcells_GOresults_df <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", header = TRUE)
gdTcells_GOresults_df$Description <- trimws(as.character(gdTcells_GOresults_df$Description))



print("CD2- gd T cells actin cytoskeleton organization Plot")
actin_cytoskeleton_organization_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "actin cytoskeleton organization") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOaco <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% actin_cytoskeleton_organization_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOaco, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOaco_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOaco_DEGs<-gdTcells_GOaco$gene
gdTcells_GOaco_DEGs <- unique(gdTcells_GOaco_DEGs)
print("length(gdTcells_GOaco_DEGs)")
length(gdTcells_GOaco_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOaco$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOaco$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOaco_breaks <- c(
    (min(gdTcells_GOaco$avg_log2FC) + (min(gdTcells_GOaco$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOaco$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOaco$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOaco$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOaco$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOaco$avg_log2FC) + (max(gdTcells_GOaco$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOaco_breaks<-round(gdTcells_GOaco_breaks, 2)  # Round normally for the rest
gdTcells_GOaco_breaks<- c(min_value, gdTcells_GOaco_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOaco_breaks")
gdTcells_GOaco_breaks
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
gdTcells_GOaco_breaks <- c(
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
gdTcells_GOaco_breaks<-round(gdTcells_GOaco_breaks, 2)  # Round normally for the rest
print("gdTcells_GOaco_breaks")
gdTcells_GOaco_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOactin_cytoskeleton_organization_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOaco, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOaco_breaks,
      labels =  gdTcells_GOaco_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'actin cytoskeleton organization' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells apoptotic process Plot")
apoptotic_process_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "apoptotic process") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOap <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% apoptotic_process_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOap, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOap_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOap_DEGs<-gdTcells_GOap$gene
gdTcells_GOap_DEGs <- unique(gdTcells_GOap_DEGs)
print("length(gdTcells_GOap_DEGs)")
length(gdTcells_GOap_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOap$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOap$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOap_breaks <- c(
    (min(gdTcells_GOap$avg_log2FC) + (min(gdTcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOap$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOap$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOap$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOap$avg_log2FC) + (max(gdTcells_GOap$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOap_breaks<-round(gdTcells_GOap_breaks, 2)  # Round normally for the rest
gdTcells_GOap_breaks<- c(min_value, gdTcells_GOap_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOap_breaks")
gdTcells_GOap_breaks
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
gdTcells_GOap_breaks <- c(
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
gdTcells_GOap_breaks<-round(gdTcells_GOap_breaks, 2)  # Round normally for the rest
print("gdTcells_GOap_breaks")
gdTcells_GOap_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOapoptotic_process_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOap, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOap_breaks,
      labels =  gdTcells_GOap_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'apoptotic process' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells apoptotic process paper Plot")
apoptotic_process_gdTcells_paper <- c("UBB","TMBIM6","SRGN","SON","SLC25A6","RPL26","RACK1","PLAC8","PARK7","MIF","GATA3","FLNA","ENO1","EMP3","CORO1A","CFL1","ANXA1","ACTB")
  gdTcells_GOap_paper <- gdTcells_day %>%
  dplyr::filter(gene %in% apoptotic_process_gdTcells_paper)
write.table(gdTcells_GOap_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOap_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOap_DEGs<-gdTcells_GOap_paper$gene
gdTcells_GOap_DEGs <- unique(gdTcells_GOap_DEGs)
print("length(gdTcells_GOap_DEGs)")
length(gdTcells_GOap_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOap_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOap_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOap_paper_breaks <- c(
    (min(gdTcells_GOap_paper$avg_log2FC) + (min(gdTcells_GOap_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOap_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOap_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOap_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOap_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOap_paper$avg_log2FC) + (max(gdTcells_GOap_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOap_paper_breaks<-round(gdTcells_GOap_paper_breaks, 2)  # Round normally for the rest
gdTcells_GOap_paper_breaks<- c(min_value, gdTcells_GOap_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOap_paper_breaks")
gdTcells_GOap_paper_breaks
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
gdTcells_GOap_paper_breaks <- c(
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
gdTcells_GOap_paper_breaks<-round(gdTcells_GOap_paper_breaks, 2)  # Round normally for the rest
print("gdTcells_GOap_paper_breaks")
gdTcells_GOap_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOapoptotic_process_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOap_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOap_paper_breaks,
      labels =  gdTcells_GOap_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'apoptotic process' in CD2- gd T cells for paper")  # Add title
dev.off()

print("CD2- gd T cells cell adhesion molecule binding Plot")
cell_adhesion_molecule_binding_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "cell adhesion molecule binding") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOcamb <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% cell_adhesion_molecule_binding_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOcamb, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOcamb_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOcamb_DEGs<-gdTcells_GOcamb$gene
gdTcells_GOcamb_DEGs <- unique(gdTcells_GOcamb_DEGs)
print("length(gdTcells_GOcamb_DEGs)")
length(gdTcells_GOcamb_DEGs)
print("gdTcells_GOcamb_DEGs")
gdTcells_GOcamb_DEGs
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOcamb$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOcamb$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOcamb_breaks <- c(
    (min(gdTcells_GOcamb$avg_log2FC) + (min(gdTcells_GOcamb$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOcamb$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOcamb$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOcamb$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOcamb$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOcamb$avg_log2FC) + (max(gdTcells_GOcamb$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOcamb_breaks<-round(gdTcells_GOcamb_breaks, 2)  # Round normally for the rest
gdTcells_GOcamb_breaks<- c(min_value, gdTcells_GOcamb_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOcamb_breaks")
gdTcells_GOcamb_breaks
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
gdTcells_GOcamb_breaks <- c(
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
gdTcells_GOcamb_breaks<-round(gdTcells_GOcamb_breaks, 2)  # Round normally for the rest
print("gdTcells_GOcamb_breaks")
gdTcells_GOcamb_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOcell_adhesion_molecule_binding_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOcamb, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOcamb_breaks,
      labels =  gdTcells_GOcamb_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cell adhesion molecule binding' in CD2- gd T cells")  # Add title
dev.off()

cell_adhesion_molecule_binding_paper <- c("ANXA1","ENO1","CAPG","HSPA8","FLNA","S100A11","PARK7", "ITGB2" ,"STAT1","SERBP1","LGALS8","CD151","RACK1")
print("CD2- gd T cells cell adhesion molecule binding paper Plot")
cell_adhesion_molecule_binding_paper <- c("ANXA1","ENO1","CAPG","HSPA8","FLNA","S100A11","PARK7", "ITGB2" ,"STAT1","SERBP1","LGALS8","CD151","RACK1")
  gdTcells_GOcamb_paper <- gdTcells_day %>%
  dplyr::filter(gene %in% cell_adhesion_molecule_binding_paper)
write.table(gdTcells_GOcamb_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOcamb_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOcamb_DEGs<-gdTcells_GOcamb_paper$gene
gdTcells_GOcamb_DEGs <- unique(gdTcells_GOcamb_DEGs)
print("length(gdTcells_GOcamb_DEGs)")
length(gdTcells_GOcamb_DEGs)
print("gdTcells_GOcamb_DEGs")
gdTcells_GOcamb_DEGs
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOcamb_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOcamb_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOcamb_paper_breaks <- c(
    (min(gdTcells_GOcamb_paper$avg_log2FC) + (min(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOcamb_paper$avg_log2FC) + (max(gdTcells_GOcamb_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOcamb_paper_breaks<-round(gdTcells_GOcamb_paper_breaks, 2)  # Round normally for the rest
gdTcells_GOcamb_paper_breaks<- c(min_value, gdTcells_GOcamb_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOcamb_paper_breaks")
gdTcells_GOcamb_paper_breaks
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
gdTcells_GOcamb_paper_breaks <- c(
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
gdTcells_GOcamb_paper_breaks<-round(gdTcells_GOcamb_paper_breaks, 2)  # Round normally for the rest
print("gdTcells_GOcamb_paper_breaks")
gdTcells_GOcamb_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOcell_adhesion_molecule_binding_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOcamb_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOcamb_paper_breaks,
      labels =  gdTcells_GOcamb_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cell adhesion molecule binding' in CD2- gd T cells for paper")  # Add title
dev.off()

print("CD2- gd T cells cellular response to endogenous stimulus Plot")
cellular_response_to_endogenous_stimulus_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "cellular response to endogenous stimulus") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOcrtes <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% cellular_response_to_endogenous_stimulus_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOcrtes, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOcrtes_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOcrtes_DEGs<-gdTcells_GOcrtes$gene
gdTcells_GOcrtes_DEGs <- unique(gdTcells_GOcrtes_DEGs)
print("length(gdTcells_GOcrtes_DEGs)")
length(gdTcells_GOcrtes_DEGs)
print("gdTcells_GOcrtes_DEGs")
gdTcells_GOcrtes_DEGs
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOcrtes$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOcrtes$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOcrtes_breaks <- c(
    (min(gdTcells_GOcrtes$avg_log2FC) + (min(gdTcells_GOcrtes$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOcrtes$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOcrtes$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOcrtes$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOcrtes$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOcrtes$avg_log2FC) + (max(gdTcells_GOcrtes$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOcrtes_breaks<-round(gdTcells_GOcrtes_breaks, 2)  # Round normally for the rest
gdTcells_GOcrtes_breaks<- c(min_value, gdTcells_GOcrtes_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOcrtes_breaks")
gdTcells_GOcrtes_breaks
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
gdTcells_GOcrtes_breaks <- c(
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
gdTcells_GOcrtes_breaks<-round(gdTcells_GOcrtes_breaks, 2)  # Round normally for the rest
print("gdTcells_GOcrtes_breaks")
gdTcells_GOcrtes_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOcellular_response_to_endogenous_stimulus_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOcrtes, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOcrtes_breaks,
      labels =  gdTcells_GOcrtes_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cellular response to endogenous stimulus' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells cellular response to endogenous stimulus paper Plot")
cellular_response_to_endogenous_stimulus_paper <- c("ACTB","ANXA1","CORO1A","CD63","HSPA8","FLNA","ARPC3","PARK7","STAT1","PPP3CA","HDAC9","PDE3B","ND3","GATA3","KLF2","RACK1","ITGA4")
#papers that prove; CD63 https://academic.oup.com/jimmunol/article/173/10/6000/8059942?login=true , HSPA8 https://jlb.onlinelibrary.wiley.com/doi/full/10.1002/JLB.3AB0420-282R , ARPC3 https://www.nature.com/articles/s41598-017-08357-4 
  gdTcells_GOcrtes_paper <- gdTcells_day %>%
  dplyr::filter(gene %in% cellular_response_to_endogenous_stimulus_paper)
write.table(gdTcells_GOcrtes_paper, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOcrtes_paper_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOcrtes_DEGs<-gdTcells_GOcrtes_paper$gene
gdTcells_GOcrtes_DEGs <- unique(gdTcells_GOcrtes_DEGs)
print("length(gdTcells_GOcrtes_DEGs)")
length(gdTcells_GOcrtes_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOcrtes_paper$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOcrtes_paper$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOcrtes_paper_breaks <- c(
    (min(gdTcells_GOcrtes_paper$avg_log2FC) + (min(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOcrtes_paper$avg_log2FC) + (max(gdTcells_GOcrtes_paper$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOcrtes_paper_breaks<-round(gdTcells_GOcrtes_paper_breaks, 2)  # Round normally for the rest
gdTcells_GOcrtes_paper_breaks<- c(min_value, gdTcells_GOcrtes_paper_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOcrtes_paper_breaks")
gdTcells_GOcrtes_paper_breaks
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
gdTcells_GOcrtes_paper_breaks <- c(
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
gdTcells_GOcrtes_paper_breaks<-round(gdTcells_GOcrtes_paper_breaks, 2)  # Round normally for the rest
print("gdTcells_GOcrtes_paper_breaks")
gdTcells_GOcrtes_paper_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOcellular_response_to_endogenous_stimulus_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOcrtes_paper, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOcrtes_paper_breaks,
      labels =  gdTcells_GOcrtes_paper_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'cellular response to endogenous stimulus' in CD2- gd T cells for paper")  # Add title
dev.off()

print("CD2- gd T cells leukocyte chemotaxis Plot")
leukocyte_chemotaxis_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "leukocyte chemotaxis") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOlc <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% leukocyte_chemotaxis_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOlc, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOlc_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOlc_DEGs<-gdTcells_GOlc$gene
gdTcells_GOlc_DEGs <- unique(gdTcells_GOlc_DEGs)
print("length(gdTcells_GOlc_DEGs)")
length(gdTcells_GOlc_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOlc$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOlc$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOlc_breaks <- c(
    (min(gdTcells_GOlc$avg_log2FC) + (min(gdTcells_GOlc$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOlc$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOlc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOlc$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOlc$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOlc$avg_log2FC) + (max(gdTcells_GOlc$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOlc_breaks<-round(gdTcells_GOlc_breaks, 2)  # Round normally for the rest
gdTcells_GOlc_breaks<- c(min_value, gdTcells_GOlc_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOlc_breaks")
gdTcells_GOlc_breaks
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
gdTcells_GOlc_breaks <- c(
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
gdTcells_GOlc_breaks<-round(gdTcells_GOlc_breaks, 2)  # Round normally for the rest
print("gdTcells_GOlc_breaks")
gdTcells_GOlc_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOleukocyte_chemotaxis_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOlc, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOlc_breaks,
      labels =  gdTcells_GOlc_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'leukocyte chemotaxis' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells lysosome Plot")
lysosome_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "lysosome") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOl <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% lysosome_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOl, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOl_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOl_DEGs<-gdTcells_GOl$gene
gdTcells_GOl_DEGs <- unique(gdTcells_GOl_DEGs)
print("length(gdTcells_GOl_DEGs)")
length(gdTcells_GOl_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOl$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOl$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOl_breaks <- c(
    (min(gdTcells_GOl$avg_log2FC) + (min(gdTcells_GOl$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOl$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOl$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOl$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOl$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOl$avg_log2FC) + (max(gdTcells_GOl$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOl_breaks<-round(gdTcells_GOl_breaks, 2)  # Round normally for the rest
gdTcells_GOl_breaks<- c(min_value, gdTcells_GOl_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOl_breaks")
gdTcells_GOl_breaks
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
gdTcells_GOl_breaks <- c(
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
gdTcells_GOl_breaks<-round(gdTcells_GOl_breaks, 2)  # Round normally for the rest
print("gdTcells_GOl_breaks")
gdTcells_GOl_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOlysosome_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOl, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOl_breaks,
      labels =  gdTcells_GOl_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'lysosome' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells lytic vacuole Plot")
lytic_vacuole_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "lytic vacuole") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOlv <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% lytic_vacuole_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOlv, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOlv_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOlv_DEGs<-gdTcells_GOlv$gene
gdTcells_GOlv_DEGs <- unique(gdTcells_GOlv_DEGs)
print("length(gdTcells_GOlv_DEGs)")
length(gdTcells_GOlv_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOlv$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOlv$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOlv_breaks <- c(
    (min(gdTcells_GOlv$avg_log2FC) + (min(gdTcells_GOlv$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOlv$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOlv$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOlv$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOlv$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOlv$avg_log2FC) + (max(gdTcells_GOlv$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOlv_breaks<-round(gdTcells_GOlv_breaks, 2)  # Round normally for the rest
gdTcells_GOlv_breaks<- c(min_value, gdTcells_GOlv_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOlv_breaks")
gdTcells_GOlv_breaks
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
gdTcells_GOlv_breaks <- c(
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
gdTcells_GOlv_breaks<-round(gdTcells_GOlv_breaks, 2)  # Round normally for the rest
print("gdTcells_GOlv_breaks")
gdTcells_GOlv_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOlytic_vacuole_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOlv, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOlv_breaks,
      labels =  gdTcells_GOlv_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'lytic vacuole' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells regulation of cell development Plot")
regulation_of_cell_development_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "regulation of cell development") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOrocd <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% regulation_of_cell_development_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOrocd, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOrocd_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOrocd_DEGs<-gdTcells_GOrocd$gene
gdTcells_GOrocd_DEGs <- unique(gdTcells_GOrocd_DEGs)
print("length(gdTcells_GOrocd_DEGs)")
length(gdTcells_GOrocd_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOrocd$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOrocd$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOrocd_breaks <- c(
    (min(gdTcells_GOrocd$avg_log2FC) + (min(gdTcells_GOrocd$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOrocd$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOrocd$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOrocd$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOrocd$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOrocd$avg_log2FC) + (max(gdTcells_GOrocd$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOrocd_breaks<-round(gdTcells_GOrocd_breaks, 2)  # Round normally for the rest
gdTcells_GOrocd_breaks<- c(min_value, gdTcells_GOrocd_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOrocd_breaks")
gdTcells_GOrocd_breaks
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
gdTcells_GOrocd_breaks <- c(
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
gdTcells_GOrocd_breaks<-round(gdTcells_GOrocd_breaks, 2)  # Round normally for the rest
print("gdTcells_GOrocd_breaks")
gdTcells_GOrocd_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOregulation_of_cell_development_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOrocd, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOrocd_breaks,
      labels =  gdTcells_GOrocd_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'regulation of cell development' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells supramolecular fiber organization Plot")
supramolecular_fiber_organization_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "supramolecular fiber organization") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOsfo <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% supramolecular_fiber_organization_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOsfo, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOsfo_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOsfo_DEGs<-gdTcells_GOsfo$gene
gdTcells_GOsfo_DEGs <- unique(gdTcells_GOsfo_DEGs)
print("length(gdTcells_GOsfo_DEGs)")
length(gdTcells_GOsfo_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOsfo$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOsfo$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOsfo_breaks <- c(
    (min(gdTcells_GOsfo$avg_log2FC) + (min(gdTcells_GOsfo$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOsfo$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOsfo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOsfo$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOsfo$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOsfo$avg_log2FC) + (max(gdTcells_GOsfo$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOsfo_breaks<-round(gdTcells_GOsfo_breaks, 2)  # Round normally for the rest
gdTcells_GOsfo_breaks<- c(min_value, gdTcells_GOsfo_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOsfo_breaks")
gdTcells_GOsfo_breaks
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
gdTcells_GOsfo_breaks <- c(
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
gdTcells_GOsfo_breaks<-round(gdTcells_GOsfo_breaks, 2)  # Round normally for the rest
print("gdTcells_GOsfo_breaks")
gdTcells_GOsfo_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOsupramolecular_fiber_organization_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOsfo, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOsfo_breaks,
      labels =  gdTcells_GOsfo_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'supramolecular fiber organization' in CD2- gd T cells")  # Add title
dev.off()
print("CD2- gd T cells actin binding Plot")
actin_binding_gdTcells <- gdTcells_GOresults_df %>% 
  dplyr::filter(Description == "actin binding") %>% 
  pull(geneID) %>% 
  strsplit("/") %>% 
  unlist()
  gdTcells_GOab <- gdTcells_day %>%
  dplyr::group_by(Human.gene.stable.ID) %>%  # Group by gene
  dplyr::filter(any(Human.gene.stable.ID %in% actin_binding_gdTcells)) %>%  # Check if any row satisfies the condition
      ungroup() 
write.table(gdTcells_GOab, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/gdTcells_GOab_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#make vector of genes
gdTcells_GOab_DEGs<-gdTcells_GOab$gene
gdTcells_GOab_DEGs <- unique(gdTcells_GOab_DEGs)
print("length(gdTcells_GOab_DEGs)")
length(gdTcells_GOab_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_GOab$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_GOab$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_GOab_breaks <- c(
    (min(gdTcells_GOab$avg_log2FC) + (min(gdTcells_GOab$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_GOab$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_GOab$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_GOab$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_GOab$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_GOab$avg_log2FC) + (max(gdTcells_GOab$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_GOab_breaks<-round(gdTcells_GOab_breaks, 2)  # Round normally for the rest
gdTcells_GOab_breaks<- c(min_value, gdTcells_GOab_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_GOab_breaks")
gdTcells_GOab_breaks
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
gdTcells_GOab_breaks <- c(
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
gdTcells_GOab_breaks<-round(gdTcells_GOab_breaks, 2)  # Round normally for the rest
print("gdTcells_GOab_breaks")
gdTcells_GOab_breaks
}

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DE_GOactin_binding_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_GOab, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_GOab_breaks,
      labels =  gdTcells_GOab_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in GO Term 'actin binding' in CD2- gd T cells")  # Add title
dev.off()

print("CD2- gd T cells DE plot")
gdTcells_DEGs <-c("ZFP36L2","UBB","SRGN","RACK1","RAC2","PPP3CA","PLAC8","GATA3","FYN","FLNA","ETS1","CORO1A","CD151","B2M","ACTB")
gdTcells_DE_plot <- gdTcells_day %>%
  dplyr::filter(gene %in% gdTcells_DEGs)
write.table(gdTcells_DE_plot, "/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/CD2neg_GD_Tcells_GO_DEGs_paper_dotPlot.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#make vector of genes
gdTcells_DE_plot_DEGs<-gdTcells_DE_plot$gene
gdTcells_DE_plot_DEGs <- unique(gdTcells_DE_plot_DEGs)
print("length(gdTcells_DE_plot_DEGs)")
length(gdTcells_DE_plot_DEGs)
#make breaks
# Pre-calculate the minimum and maximum values
min_value <- min(gdTcells_DE_plot$avg_log2FC, na.rm = TRUE)
#round the min_value down to 2 decimal places
min_value <- ceiling(min_value * 100) / 100
max_value <- max(gdTcells_DE_plot$avg_log2FC, na.rm = TRUE)
#round the max_value down to 2 decimal places
max_value <- floor(max_value * 100) / 100
#if min_value is <= 0, them make breaks this way
if (min_value <= 0) {
gdTcells_DE_plot_breaks <- c(
    (min(gdTcells_DE_plot$avg_log2FC) + (min(gdTcells_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between negative  midpoint and min
    (min(gdTcells_DE_plot$avg_log2FC) + 0) / 2,  # Midpoint between min and 0
    (0 + (min(gdTcells_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and negative midpoint
    (0 + (max(gdTcells_DE_plot$avg_log2FC) + 0) / 2) / 2,  # Break between 0 and positive midpoint
    (max(gdTcells_DE_plot$avg_log2FC) + 0) / 2,  # Positive midpoint
    (max(gdTcells_DE_plot$avg_log2FC) + (max(gdTcells_DE_plot$avg_log2FC) + 0) / 2) / 2 # Break between positive midpoint and max
    )
gdTcells_DE_plot_breaks<-round(gdTcells_DE_plot_breaks, 2)  # Round normally for the rest
gdTcells_DE_plot_breaks<- c(min_value, gdTcells_DE_plot_breaks, max_value, 0)  # Add min and max to the breaks    
print("gdTcells_DE_plot_breaks")
gdTcells_DE_plot_breaks
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
gdTcells_DE_plot_breaks <- c(
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
gdTcells_DE_plot_breaks<-round(gdTcells_DE_plot_breaks, 2)  # Round normally for the rest
print("gdTcells_DE_plot_breaks")
gdTcells_DE_plot_breaks
} 

pdf("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/GO_DEG_Plots/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_gdTcells_DEGs_paper_dotPlotCSP.pdf", width=5 ,height=10)
ggplot( gdTcells_DE_plot, aes(x = Comparison, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21) +
scale_fill_gradient2(
    low = "blue4",       # Color for the lowest values
    mid = "white",      # Color for the midpoint (0)
    high = "red4",       # Color for the highest values
    midpoint = 0,       
    breaks =  gdTcells_DE_plot_breaks,
      labels =  gdTcells_DE_plot_breaks, 
      guide = guide_colorbar(reverse = FALSE)              # Keep the legend order normal
  ) +   scale_size(range = c(5,15))+  # Round labels to the nearest whole number

 scale_x_discrete(breaks = all_Comparisons ) +  # Ensure all clusters are labeled
scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +  # Wrap text to 50 characters per line
theme_classic() +  # Use classic theme
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 315, hjust = 0),  # Increase size and set color to black
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
    labs(x = "Comparison", y = "Gene", fill = expression("Avg Log"[2]*"FC"), size = "Percent Expressing") + 
    ggtitle( "DEGs in CD2- gd T cells")  # Add title
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