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


print("Go on Monocytes")
#Read in list of all genes Monocytes for GO analysis
Monocytes_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_bkg.txt", sep = "\t", header = TRUE)
Monocytes_bkg<-Monocytes_bkg$Human.gene.stable.ID
#read in C0_1v0_0 significant genes for Monocytes
C0_1v0_0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
C0_1v0_0<-distinct(C0_1v0_0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(C0_1v0_0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(C0_1v0_0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

C0_1v0_0_2 <- C0_1v0_0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
C0_1v0_0_2_list<-C0_1v0_0$avg_log2FC
# name the vector
names(C0_1v0_0_2_list) <-C0_1v0_0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
C0_1v0_0_2_list = sort(C0_1v0_0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
C0_1v0_0_2_list_pos<-C0_1v0_0_2_list[C0_1v0_0_2_list > 0]
C0_1v0_0_2_list_neg<-C0_1v0_0_2_list[C0_1v0_0_2_list < 0]
#Postive genes
C0_1v0_0_pos_enrichGO <- enrichGO(gene=names(C0_1v0_0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(C0_1v0_0_pos_enrichGO) && nrow(C0_1v0_0_pos_enrichGO) > 0) {
C0_1v0_0_pos_enrichGO_s<-saveRDS(C0_1v0_0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_2026_04_14.rds")
C0_1v0_0_pos_enrichGO_df<-as.data.frame(C0_1v0_0_pos_enrichGO)
C0_1v0_0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(C0_1v0_0_pos_enrichGO_df$GeneRatio))
C0_1v0_0_pos_enrichGO_df$Celltype <- "Monocytes"
C0_1v0_0_pos_enrichGO_df$Comparsion <- "0_1 vs 0_0"
colnames(C0_1v0_0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C0_1v0_0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C0_1v0_0_pos_enrichGO_simp<-clusterProfiler::simplify(C0_1v0_0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}
if (!is.null(C0_1v0_0_pos_enrichGO_simp) && nrow(C0_1v0_0_pos_enrichGO_simp) > 0) {
saveRDS(C0_1v0_0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_simplified_2026_04_14.rds")
C0_1v0_0_pos_enrichGO_simp_df<-as.data.frame(C0_1v0_0_pos_enrichGO_simp)
C0_1v0_0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C0_1v0_0_pos_enrichGO_simp_df$GeneRatio))
colnames(C0_1v0_0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C0_1v0_0_pos_enrichGO_simp_df$Celltype <- "Monocytes"
C0_1v0_0_pos_enrichGO_simp_df$Comparsion <- "0_1 vs 0_0"
write.table(C0_1v0_0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
C0_1v0_0_neg_enrichGO <- enrichGO(gene=names(C0_1v0_0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(C0_1v0_0_neg_enrichGO) && nrow(C0_1v0_0_neg_enrichGO) > 0) {
## Save the results
C0_1v0_0_neg_enrichGO_s<-saveRDS(C0_1v0_0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_neg_enrichGO_2026_04_14.rds")
C0_1v0_0_neg_enrichGO_df<-as.data.frame(C0_1v0_0_neg_enrichGO)
C0_1v0_0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(C0_1v0_0_neg_enrichGO_df$GeneRatio))
C0_1v0_0_neg_enrichGO_df$Celltype <- "Monocytes"
C0_1v0_0_neg_enrichGO_df$Comparsion <- "0_1 vs 0_0"
colnames(C0_1v0_0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C0_1v0_0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C0_1v0_0_neg_enrichGO_simp<-clusterProfiler::simplify(C0_1v0_0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}
if (!is.null(C0_1v0_0_neg_enrichGO_simp) && nrow(C0_1v0_0_neg_enrichGO_simp) > 0) {
saveRDS(C0_1v0_0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_neg_enrichGO_simplified_2026_04_14.rds")
C0_1v0_0_neg_enrichGO_simp_df<-as.data.frame(C0_1v0_0_neg_enrichGO_simp)
C0_1v0_0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C0_1v0_0_neg_enrichGO_simp_df$GeneRatio))
colnames(C0_1v0_0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C0_1v0_0_neg_enrichGO_simp_df$Celltype <- "Monocytes"
C0_1v0_0_neg_enrichGO_simp_df$Comparsion <- "0_1 vs 0_0"
write.table(C0_1v0_0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

##run D8 vs D2
#read in C0v1 significant genes for Monocytes
C0v1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
C0v1<-distinct(C0v1)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(C0v1$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(C0v1$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

C0v1_2 <- C0v1 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
C0v1_2_list<-C0v1$avg_log2FC
# name the vector
names(C0v1_2_list) <-C0v1_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
C0v1_2_list = sort(C0v1_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
C0v1_2_list_pos<-C0v1_2_list[C0v1_2_list > 0]
C0v1_2_list_neg<-C0v1_2_list[C0v1_2_list < 0]
#Postive genes
C0v1_pos_enrichGO <- enrichGO(gene=names(C0v1_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(C0v1_pos_enrichGO) && nrow(C0v1_pos_enrichGO) > 0) {
C0v1_pos_enrichGO_s<-saveRDS(C0v1_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_2026_04_14.rds")
C0v1_pos_enrichGO_df<-as.data.frame(C0v1_pos_enrichGO)
C0v1_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(C0v1_pos_enrichGO_df$GeneRatio))
C0v1_pos_enrichGO_df$Celltype <- "Monocytes"
C0v1_pos_enrichGO_df$Comparsion <- "0_0 vs 0_1"
colnames(C0v1_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C0v1_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C0v1_pos_enrichGO_simp<-clusterProfiler::simplify(C0v1_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C0v1_pos_enrichGO_simp) && nrow(C0v1_pos_enrichGO_simp) > 0) {
saveRDS(C0v1_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_simplified_2026_04_14.rds")
C0v1_pos_enrichGO_simp_df<-as.data.frame(C0v1_pos_enrichGO_simp)
C0v1_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C0v1_pos_enrichGO_simp_df$GeneRatio))
colnames(C0v1_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C0v1_pos_enrichGO_simp_df$Celltype <- "Monocytes"
C0v1_pos_enrichGO_simp_df$Comparsion <- "0_0 vs 0_1"
write.table(C0v1_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
C0v1_neg_enrichGO <- enrichGO(gene=names(C0v1_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(C0v1_neg_enrichGO) && nrow(C0v1_neg_enrichGO) > 0) {
## Save the results
C0v1_neg_enrichGO_s<-saveRDS(C0v1_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_neg_enrichGO_2026_04_14.rds")
C0v1_neg_enrichGO_df<-as.data.frame(C0v1_neg_enrichGO)
C0v1_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(C0v1_neg_enrichGO_df$GeneRatio))
C0v1_neg_enrichGO_df$Celltype <- "Monocytes"
C0v1_neg_enrichGO_df$Comparsion <- "0_0 vs 0_1"
colnames(C0v1_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C0v1_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C0v1_neg_enrichGO_simp<-clusterProfiler::simplify(C0v1_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C0v1_neg_enrichGO_simp) && nrow(C0v1_neg_enrichGO_simp) > 0) {
saveRDS(C0v1_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_neg_enrichGO_simplified_2026_04_14.rds")
C0v1_neg_enrichGO_simp_df<-as.data.frame(C0v1_neg_enrichGO_simp)
C0v1_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C0v1_neg_enrichGO_simp_df$GeneRatio))
colnames(C0v1_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C0v1_neg_enrichGO_simp_df$Celltype <- "Monocytes"
C0v1_neg_enrichGO_simp_df$Comparsion <- "0_0 vs 0_1"
write.table(C0v1_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

print("Go on NK cells")
#Read in list of all genes NKcells for GO analysis
NKcells_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_bkg.txt", sep = "\t", header = TRUE)
NKcells_bkg<-NKcells_bkg$Human.gene.stable.ID
#read in C4_1v4_0 significant genes for NKcells
C4_1v4_0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
C4_1v4_0<-distinct(C4_1v4_0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(C4_1v4_0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(C4_1v4_0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

C4_1v4_0_2 <- C4_1v4_0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
C4_1v4_0_2_list<-C4_1v4_0$avg_log2FC
# name the vector
names(C4_1v4_0_2_list) <-C4_1v4_0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
C4_1v4_0_2_list = sort(C4_1v4_0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
C4_1v4_0_2_list_pos<-C4_1v4_0_2_list[C4_1v4_0_2_list > 0]
C4_1v4_0_2_list_neg<-C4_1v4_0_2_list[C4_1v4_0_2_list < 0]
#Postive genes
C4_1v4_0_pos_enrichGO <- enrichGO(gene=names(C4_1v4_0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(C4_1v4_0_pos_enrichGO) && nrow(C4_1v4_0_pos_enrichGO) > 0) {
C4_1v4_0_pos_enrichGO_s<-saveRDS(C4_1v4_0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_2026_04_14.rds")
C4_1v4_0_pos_enrichGO_df<-as.data.frame(C4_1v4_0_pos_enrichGO)
C4_1v4_0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(C4_1v4_0_pos_enrichGO_df$GeneRatio))
C4_1v4_0_pos_enrichGO_df$Celltype <- "NK cells"
C4_1v4_0_pos_enrichGO_df$Comparsion <- "4_1 vs 4_0"
colnames(C4_1v4_0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C4_1v4_0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_1v4_0_pos_enrichGO_simp<-clusterProfiler::simplify(C4_1v4_0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C4_1v4_0_pos_enrichGO_simp) && nrow(C4_1v4_0_pos_enrichGO_simp) > 0) {
saveRDS(C4_1v4_0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_simplified_2026_04_14.rds")
C4_1v4_0_pos_enrichGO_simp_df<-as.data.frame(C4_1v4_0_pos_enrichGO_simp)
C4_1v4_0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C4_1v4_0_pos_enrichGO_simp_df$GeneRatio))
colnames(C4_1v4_0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C4_1v4_0_pos_enrichGO_simp_df$Celltype <- "NK cells"
C4_1v4_0_pos_enrichGO_simp_df$Comparsion <- "4_1 vs 4_0"
write.table(C4_1v4_0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
C4_1v4_0_neg_enrichGO <- enrichGO(gene=names(C4_1v4_0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(C4_1v4_0_neg_enrichGO) && nrow(C4_1v4_0_neg_enrichGO) > 0) {
## Save the results
C4_1v4_0_neg_enrichGO_s<-saveRDS(C4_1v4_0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_neg_enrichGO_2026_04_14.rds")
C4_1v4_0_neg_enrichGO_df<-as.data.frame(C4_1v4_0_neg_enrichGO)
C4_1v4_0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(C4_1v4_0_neg_enrichGO_df$GeneRatio))
C4_1v4_0_neg_enrichGO_df$Celltype <- "NK cells"
C4_1v4_0_neg_enrichGO_df$Comparsion <- "4_1 vs 4_0"
colnames(C4_1v4_0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C4_1v4_0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_1v4_0_neg_enrichGO_simp<-clusterProfiler::simplify(C4_1v4_0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C4_1v4_0_neg_enrichGO_simp) && nrow(C4_1v4_0_neg_enrichGO_simp) > 0) {
saveRDS(C4_1v4_0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_neg_enrichGO_simplified_2026_04_14.rds")
C4_1v4_0_neg_enrichGO_simp_df<-as.data.frame(C4_1v4_0_neg_enrichGO_simp)
C4_1v4_0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C4_1v4_0_neg_enrichGO_simp_df$GeneRatio))
colnames(C4_1v4_0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C4_1v4_0_neg_enrichGO_simp_df$Celltype <- "NK cells"
C4_1v4_0_neg_enrichGO_simp_df$Comparsion <- "4_1 vs 4_0"
write.table(C4_1v4_0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

##run D8 vs D2
#read in C0v1 significant genes for NKcells
C4_0v4_1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
C4_0v4_1<-distinct(C4_0v4_1)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(C4_0v4_1$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(C4_0v4_1$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

C4_0v4_1 <- C4_0v4_1 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
C4_0v4_1_list<-C4_0v4_1$avg_log2FC
# name the vector
names(C4_0v4_1_list) <- C4_0v4_1$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
C4_0v4_1_list = sort(C4_0v4_1_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
C4_0v4_1_list_pos<-C4_0v4_1_list[C4_0v4_1_list > 0]
C4_0v4_1_list_neg<-C4_0v4_1_list[C4_0v4_1_list < 0]   
#Postive genes
C4_0v4_1_pos_enrichGO <- enrichGO(gene=names(C4_0v4_1_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(C4_0v4_1_pos_enrichGO) && nrow(C4_0v4_1_pos_enrichGO) > 0) {
C4_0v4_1_pos_enrichGO_s<-saveRDS(C4_0v4_1_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_2026_04_14.rds")
C4_0v4_1_pos_enrichGO_df<-as.data.frame(C4_0v4_1_pos_enrichGO)
C4_0v4_1_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(C4_0v4_1_pos_enrichGO_df$GeneRatio))
C4_0v4_1_pos_enrichGO_df$Celltype <- "NK cells"
C4_0v4_1_pos_enrichGO_df$Comparsion <- "4_0 vs 4_1"
colnames(C4_0v4_1_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C4_0v4_1_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_0v4_1_pos_enrichGO_simp<-clusterProfiler::simplify(C4_0v4_1_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C4_0v4_1_pos_enrichGO_simp) && nrow(C4_0v4_1_pos_enrichGO_simp) > 0) {
saveRDS(C4_0v4_1_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_simplified_2026_04_14.rds")
C4_0v4_1_pos_enrichGO_simp_df<-as.data.frame(C4_0v4_1_pos_enrichGO_simp)
C4_0v4_1_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C4_0v4_1_pos_enrichGO_simp_df$GeneRatio))
colnames(C4_0v4_1_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C4_0v4_1_pos_enrichGO_simp_df$Celltype <- "NK cells"
C4_0v4_1_pos_enrichGO_simp_df$Comparsion <- "4_0 vs 4_1"
write.table(C4_0v4_1_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
C4_0v4_1_neg_enrichGO <- enrichGO(gene=names(C0v1_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(C4_0v4_1_neg_enrichGO) && nrow(C4_0v4_1_neg_enrichGO) > 0) {
## Save the results
C4_0v4_1_neg_enrichGO_s<-saveRDS(C4_0v4_1_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_neg_enrichGO_2026_04_14.rds")
C4_0v4_1_neg_enrichGO_df<-as.data.frame(C4_0v4_1_neg_enrichGO)
C4_0v4_1_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(C4_0v4_1_neg_enrichGO_df$GeneRatio))
C4_0v4_1_neg_enrichGO_df$Celltype <- "NK cells"
C4_0v4_1_neg_enrichGO_df$Comparsion <- "4_0 vs 4_1"
colnames(C4_0v4_1_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(C4_0v4_1_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_0v4_1_neg_enrichGO_simp<-clusterProfiler::simplify(C4_0v4_1_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(C4_0v4_1_neg_enrichGO_simp) && nrow(C4_0v4_1_neg_enrichGO_simp) > 0) {
saveRDS(C4_0v4_1_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_neg_enrichGO_simplified_2026_04_14.rds")
C4_0v4_1_neg_enrichGO_simp_df<-as.data.frame(C4_0v4_1_neg_enrichGO_simp)
C4_0v4_1_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(C4_0v4_1_neg_enrichGO_simp_df$GeneRatio))
colnames(C4_0v4_1_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
C4_0v4_1_neg_enrichGO_simp_df$Celltype <- "NK cells"
C4_0v4_1_neg_enrichGO_simp_df$Comparsion <- "4_0 vs 4_1"  
write.table(C4_0v4_1_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}






C0_1v0_0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
C0_0v0_01 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_all <- rbind(C0_1v0_0, C0_0v0_01)
write.table(Monocyte_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_Monocytes_subcluster_pairwise_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C0_1v0_0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_1v0_0_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
C0_0v0_01_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/Monocytes_subcluster_0_0v0_1_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_simp <- rbind(C0_1v0_0_simp, C0_0v0_01_simp)
write.table(Monocyte_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_Monocytes_subcluster_pairwise_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_1v4_0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
C4_0v4_1 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_all <- rbind(C4_1v4_0, C4_0v4_1)
write.table(NKcell_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_NKcells_subcluster_pairwise_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
C4_1v4_0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_1v4_0_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
C4_0v4_1_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/NKcells_subcluster_4_0v4_1_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_simp <- rbind(C4_1v4_0_simp, C4_0v4_1_simp)
write.table(NKcell_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_subcluster/All_NKcells_subcluster_pairwise_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)


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
# [124] crayon_1.5.3            future.apply_1.11.3     rematch2_2.1.2         
# [127] plyr_1.8.9              fs_1.6.6                ggbeeswarm_0.7.2       
# [130] stringi_1.8.7           viridisLite_0.4.2       deldir_2.0-4           
# [133] BiocParallel_1.36.0     babelgene_22.9          Biostrings_2.70.3      
# [136] lazyeval_0.2.2          spatstat.geom_3.3-5     GOSemSim_2.28.1        
# [139] RcppHNSW_0.6.0          hms_1.1.3               patchwork_1.3.0        
# [142] bit64_4.6.0-1           future_1.34.0           KEGGREST_1.42.0        
# [145] shiny_1.10.0            ROCR_1.0-11             igraph_2.1.4           
# [148] memoise_2.0.1           ggtree_3.10.1           fastmatch_1.1-6        
# [151] bit_4.6.0               gson_0.1.0              ape_5.8-1              