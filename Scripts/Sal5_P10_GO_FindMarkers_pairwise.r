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

print("Go on cDCs")
#Read in list of all genes cDCs for GO analysis
cDCs_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_bkg.txt", sep = "\t", header = TRUE)
cDCs_bkg<-cDCs_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for cDCs
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
if (!is.null(D2v0_2_list_pos) && length(D2v0_2_list_pos) > 0) {
print("D2v0_pos_enrichGO")
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = cDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
}
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "cDCs"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)  
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "cDCs"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
if (!is.null(D2v0_2_list_neg) && length(D2v0_2_list_neg) > 0) {
print("D2v0_neg_enrichGO")
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = cDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
}

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "cDCs"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}
if (exists("D2v0_neg_enrichGO_simp")) {
    if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
        saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_neg_enrichGO_simplified_2026_04_14.rds")
        D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
        D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
        colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
        D2v0_neg_enrichGO_simp_df$Celltype <- "cDCs"
        D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
        write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
    }
}


##run D8 vs D2
#read in D8v2 significant genes for cDCs
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
if (!is.null(D8v2_2_list_pos) && length(D8v2_2_list_pos) > 0) {
print("D8v2_pos_enrichGO")
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = cDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
}
## Save the results
if (exists("D8v2_pos_enrichGO")) {
    if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {k
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "cDCs"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
print("D8v2_pos_enrichGO_df")
print(head(D8v2_pos_enrichGO_df))
  # Write table with error handling
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
 }
}

if (exists("D8v2_pos_enrichGO_simp")) {
  if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {

saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "cDCs"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
}
#Negative genes
if (!is.null(D8v2_2_list_neg) && length(D8v2_2_list_neg) > 0) {
print("D8v2_neg_enrichGO")
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = cDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
}

if (exists("D8v2_neg_enrichGO") && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "cDCs"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (exists("D8v2_neg_enrichGO_simp") && !is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "cDCs"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_cDCs_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#0 DEGs found in D8 vs D0 for cDCs, so we will not run GO for that comparison
rm()


print("Go on CD4+ AB T cells")
#Read in list of all genes CD4Tcells for GO analysis
CD4Tcells_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_bkg.txt", sep = "\t", header = TRUE)
CD4Tcells_bkg<-CD4Tcells_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for CD4Tcells
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "CD4Tcells"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)") 
D2v0_pos_enrichGO_simp_df$Celltype <- "CD4Tcells"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"

write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "CD4Tcells"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio)) 
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "CD4Tcells"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
##run D8 vs D2
#read in D8v2 significant genes for CD4Tcells
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "CD4Tcells"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "CD4Tcells"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "CD4Tcells"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)  
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "CD4Tcells"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"  
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
###run D8 vs D0
#read in D8v0 significant genes for CD4Tcells
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "CD4Tcells"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "CD4Tcells"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD4Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "CD4Tcells"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "CD4Tcells"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()

print("Go on Monocytes")
#Read in list of all genes Monocytes for GO analysis
Monocytes_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_bkg.txt", sep = "\t", header = TRUE)
Monocytes_bkg<-Monocytes_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for Monocytes
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "Monocytes"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}
if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "Monocytes"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "Monocytes"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}
if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "Monocytes"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

##run D8 vs D2
#read in D8v2 significant genes for Monocytes
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "Monocytes"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "Monocytes"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "Monocytes"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "Monocytes"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

###run D8 vs D0
#read in D8v0 significant genes for Monocytes
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "Monocytes"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "Monocytes"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Monocytes_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "Monocytes"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "Monocytes"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()

print("Go on NK cells")
#Read in list of all genes NKcells for GO analysis
NKcells_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_bkg.txt", sep = "\t", header = TRUE)
NKcells_bkg<-NKcells_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for NKcells
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "NK cells"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "NK cells"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "NK cells"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "NK cells"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

##run D8 vs D2
#read in D8v2 significant genes for NKcells
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "NK cells"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "NK cells"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "NK cells"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "NK cells"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"  
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

###run D8 vs D0
#read in D8v0 significant genes for NKcells
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "NK cells"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "NK cells"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)  
}

#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = NKcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "NK cells"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "NK cells"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()

print("Go on B cells")
#Read in list of all genes Bcells for GO analysis
Bcells_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_bkg.txt", sep = "\t", header = TRUE)
Bcells_bkg<-Bcells_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for Bcells
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "B cells"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)  
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "B cells"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "B cells"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "B cells"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
##run D8 vs D2
#read in D8v2 significant genes for Bcells
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "B cells"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "B cells"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "B cells"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "B cells"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
###run D8 vs D0
#read in D8v0 significant genes for Bcells
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "B cells"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)  
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "B cells"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = Bcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "B cells"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "B cells"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()

print("Go on CD2- gd T cells")
#Read in list of all genes gCD2neg_GD_Tcells for GO analysis
CD2neg_GD_Tcells_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_bkg.txt", sep = "\t", header = TRUE)
CD2neg_GD_Tcells_bkg<-CD2neg_GD_Tcells_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for CD2neg_GD_Tcells
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "CD2- gd T cells"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "CD2- gd T cells"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
##run D8 vs D2
#read in D8v2 significant genes for CD2neg_GD_Tcells
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "CD2- gd T cells"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2" 
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "CD2- gd T cells"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)  
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
###run D8 vs D0
#read in D8v0 significant genes for CD2neg_GD_Tcells
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "CD2- gd T cells"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}

#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = CD2neg_GD_Tcells_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "CD2- gd T cells"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "CD2- gd T cells"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()


print("Go on ASCs")
#Read in list of all genes ASCs for GO analysis
ASCs_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_bkg.txt", sep = "\t", header = TRUE)
ASCs_bkg<-ASCs_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for ASCs
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "ASCs"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp <- clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)  
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "ASCs"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "ASCs"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_neg_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp <- clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"  
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_neg_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
##run D8 vs D2
#read in D8v2 significant genes for ASCs
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "ASCs"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp <- clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "ASCs"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "ASCs"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_neg_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")  
D8v2_neg_enrichGO_simp_df$Celltype <- "ASCs"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_neg_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}

###run D8 vs D0
#read in D8v0 significant genes for ASCs
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "ASCs"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)  
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "ASCs"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = ASCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "ASCs"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_neg_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "ASCs"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()

print("Go on pDCs")
#Read in list of all genes pDCs for GO analysis
pDCs_bkg<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_bkg.txt", sep = "\t", header = TRUE)
pDCs_bkg<-pDCs_bkg$Human.gene.stable.ID
#read in D2v0 significant genes for pDCs
D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D2v0<-distinct(D2v0)
#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D2v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D2v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D2v0_2 <- D2v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D2v0_2_list<-D2v0$avg_log2FC
# name the vector
names(D2v0_2_list) <-D2v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D2v0_2_list = sort(D2v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D2v0_2_list_pos<-D2v0_2_list[D2v0_2_list > 0]
D2v0_2_list_neg<-D2v0_2_list[D2v0_2_list < 0]
#Postive genes
D2v0_pos_enrichGO <- enrichGO(gene=names(D2v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D2v0_pos_enrichGO) && nrow(D2v0_pos_enrichGO) > 0) {
D2v0_pos_enrichGO_s<-saveRDS(D2v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_2026_04_14.rds")
D2v0_pos_enrichGO_df<-as.data.frame(D2v0_pos_enrichGO)
D2v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_df$GeneRatio))
D2v0_pos_enrichGO_df$Celltype <- "pDCs"
D2v0_pos_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_pos_enrichGO_simp<-clusterProfiler::simplify(D2v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_pos_enrichGO_simp) && nrow(D2v0_pos_enrichGO_simp) > 0) {
saveRDS(D2v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_simplified_2026_04_14.rds")
D2v0_pos_enrichGO_simp_df<-as.data.frame(D2v0_pos_enrichGO_simp)
D2v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D2v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_pos_enrichGO_simp_df$Celltype <- "pDCs"
D2v0_pos_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
#Negative genes
D2v0_neg_enrichGO <- enrichGO(gene=names(D2v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D2v0_neg_enrichGO) && nrow(D2v0_neg_enrichGO) > 0) {
## Save the results
D2v0_neg_enrichGO_s<-saveRDS(D2v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_neg_enrichGO_2026_04_14.rds")
D2v0_neg_enrichGO_df<-as.data.frame(D2v0_neg_enrichGO)
D2v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_df$GeneRatio))
D2v0_neg_enrichGO_df$Celltype <- "pDCs"
D2v0_neg_enrichGO_df$Comparison <- "D2 vs D0"
colnames(D2v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D2v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D2v0_neg_enrichGO_simp<-clusterProfiler::simplify(D2v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D2v0_neg_enrichGO_simp) && nrow(D2v0_neg_enrichGO_simp) > 0) {
saveRDS(D2v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_neg_enrichGO_simplified_2026_04_14.rds")
D2v0_neg_enrichGO_simp_df<-as.data.frame(D2v0_neg_enrichGO_simp)  
D2v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D2v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D2v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D2v0_neg_enrichGO_simp_df$Celltype <- "pDCs"
D2v0_neg_enrichGO_simp_df$Comparison <- "D2 vs D0"
write.table(D2v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_neg_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
##run D8 vs D2
#read in D8v2 significant genes for pDCs
D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v2<-distinct(D8v2)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v2$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v2$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v2_2 <- D8v2 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v2_2_list<-D8v2$avg_log2FC
# name the vector
names(D8v2_2_list) <-D8v2_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v2_2_list = sort(D8v2_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v2_2_list_pos<-D8v2_2_list[D8v2_2_list > 0]
D8v2_2_list_neg<-D8v2_2_list[D8v2_2_list < 0]
#Postive genes
D8v2_pos_enrichGO <- enrichGO(gene=names(D8v2_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v2_pos_enrichGO) && nrow(D8v2_pos_enrichGO) > 0) {
D8v2_pos_enrichGO_s<-saveRDS(D8v2_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_2026_04_14.rds")
D8v2_pos_enrichGO_df<-as.data.frame(D8v2_pos_enrichGO)
D8v2_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_df$GeneRatio))
D8v2_pos_enrichGO_df$Celltype <- "pDCs"
D8v2_pos_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v2_pos_enrichGO_simp<-clusterProfiler::simplify(D8v2_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_pos_enrichGO_simp) && nrow(D8v2_pos_enrichGO_simp) > 0) {
saveRDS(D8v2_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_simplified_2026_04_14.rds")
D8v2_pos_enrichGO_simp_df<-as.data.frame(D8v2_pos_enrichGO_simp)
D8v2_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v2_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_pos_enrichGO_simp_df$Celltype <- "pDCs"
D8v2_pos_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_simplified_2026_04_14.txt",sep = "\t", row.names = FALSE)
}
#Negative genes
D8v2_neg_enrichGO <- enrichGO(gene=names(D8v2_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v2_neg_enrichGO) && nrow(D8v2_neg_enrichGO) > 0) {
## Save the results
D8v2_neg_enrichGO_s<-saveRDS(D8v2_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_neg_enrichGO_2026_04_14.rds")
D8v2_neg_enrichGO_df<-as.data.frame(D8v2_neg_enrichGO)
D8v2_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_df$GeneRatio))
D8v2_neg_enrichGO_df$Celltype <- "pDCs"
D8v2_neg_enrichGO_df$Comparison <- "D8 vs D2"
colnames(D8v2_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v2_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_neg_enrichGO_2026_04_14.txt",sep = "\t", row.names = FALSE)
D8v2_neg_enrichGO_simp<-clusterProfiler::simplify(D8v2_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v2_neg_enrichGO_simp) && nrow(D8v2_neg_enrichGO_simp) > 0) {
saveRDS(D8v2_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_neg_enrichGO_simplified_2026_04_14.rds")
D8v2_neg_enrichGO_simp_df<-as.data.frame(D8v2_neg_enrichGO_simp)
D8v2_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v2_neg_enrichGO_simp_df$GeneRatio)) 
colnames(D8v2_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v2_neg_enrichGO_simp_df$Celltype <- "pDCs"
D8v2_neg_enrichGO_simp_df$Comparison <- "D8 vs D2"
write.table(D8v2_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
###run D8 vs D0
#read in D8v0 significant genes for pDCs
D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep = "\t", header = TRUE)
D8v0<-distinct(D8v0)

#See the max and min avg_log2FC, If there are negative values, we will need to split the list into positive and negative avg_log2FC values
max_avg_log2FC <- max(D8v0$avg_log2FC, na.rm = TRUE)
min_avg_log2FC <- min(D8v0$avg_log2FC, na.rm = TRUE)
print(paste("Max avg_log2FC:", max_avg_log2FC))
print(paste("Min avg_log2FC:", min_avg_log2FC))

D8v0_2 <- D8v0 %>% 
  filter(p_val_adj < 0.05 & (avg_log2FC> 0.25 | avg_log2FC < -0.25))
# we want the log2 fold change 
D8v0_2_list<-D8v0$avg_log2FC
# name the vector
names(D8v0_2_list) <-D8v0_2$Human.gene.stable.ID
# sort the list in decreasing order (required for clusterProfiler)
D8v0_2_list = sort(D8v0_2_list, decreasing = TRUE)
# Make a positive & negative list based on the logFC
D8v0_2_list_pos<-D8v0_2_list[D8v0_2_list > 0]
D8v0_2_list_neg<-D8v0_2_list[D8v0_2_list < 0]
#Postive genes
D8v0_pos_enrichGO <- enrichGO(gene=names(D8v0_2_list_pos),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  
## Save the results
if (!is.null(D8v0_pos_enrichGO) && nrow(D8v0_pos_enrichGO) > 0) {
D8v0_pos_enrichGO_s<-saveRDS(D8v0_pos_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_2026_04_14.rds")
D8v0_pos_enrichGO_df<-as.data.frame(D8v0_pos_enrichGO)
D8v0_pos_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_df$GeneRatio))
D8v0_pos_enrichGO_df$Celltype <- "pDCs"
D8v0_pos_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_pos_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_pos_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_pos_enrichGO_simp<-clusterProfiler::simplify(D8v0_pos_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_pos_enrichGO_simp) && nrow(D8v0_pos_enrichGO_simp) > 0) {
saveRDS(D8v0_pos_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_simplified_2026_04_14.rds")
D8v0_pos_enrichGO_simp_df<-as.data.frame(D8v0_pos_enrichGO_simp)
D8v0_pos_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_pos_enrichGO_simp_df$GeneRatio))
colnames(D8v0_pos_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_pos_enrichGO_simp_df$Celltype <- "pDCs"
D8v0_pos_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_pos_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
#Negative genes
D8v0_neg_enrichGO <- enrichGO(gene=names(D8v0_2_list_neg),OrgDb=org.Hs.eg.db,keyType = "ENSEMBL",ont = "ALL",pvalueCutoff = 1,pAdjustMethod = "BH",universe = pDCs_bkg,qvalueCutoff = 1,minGSSize = 10,maxGSSize = 1118964)  

if (!is.null(D8v0_neg_enrichGO) && nrow(D8v0_neg_enrichGO) > 0) {
## Save the results
D8v0_neg_enrichGO_s<-saveRDS(D8v0_neg_enrichGO, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_neg_enrichGO_2026_04_14.rds")
D8v0_neg_enrichGO_df<-as.data.frame(D8v0_neg_enrichGO)
D8v0_neg_enrichGO_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_df$GeneRatio))
D8v0_neg_enrichGO_df$Celltype <- "pDCs"
D8v0_neg_enrichGO_df$Comparison <- "D8 vs D0"
colnames(D8v0_neg_enrichGO_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
write.table(D8v0_neg_enrichGO_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_neg_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
D8v0_neg_enrichGO_simp<-clusterProfiler::simplify(D8v0_neg_enrichGO, cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)
}

if (!is.null(D8v0_neg_enrichGO_simp) && nrow(D8v0_neg_enrichGO_simp) > 0) {
saveRDS(D8v0_neg_enrichGO_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_neg_enrichGO_simplified_2026_04_14.rds")
D8v0_neg_enrichGO_simp_df<-as.data.frame(D8v0_neg_enrichGO_simp)
D8v0_neg_enrichGO_simp_df$GeneRatio<-paste0("'", as.character(D8v0_neg_enrichGO_simp_df$GeneRatio))
colnames(D8v0_neg_enrichGO_simp_df)[4]<-c("GeneRatio(' added to read in excel correctly)")
D8v0_neg_enrichGO_simp_df$Celltype <- "pDCs"
D8v0_neg_enrichGO_simp_df$Comparison <- "D8 vs D0"
write.table(D8v0_neg_enrichGO_simp_df, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_neg_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
}
rm()


Monocyte_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_all <- rbind(Monocyte_D2v0, Monocyte_D8v0, Monocyte_D8v2)
write.table(Monocyte_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Monocytes_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
Monocyte_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Monocyte_simp <- rbind(Monocyte_D2v0_simp, Monocyte_D8v0_simp, Monocyte_D8v2_simp)
write.table(Monocyte_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Monocytes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
NKcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_all <- rbind(NKcell_D2v0, NKcell_D8v0, NKcell_D8v2)
write.table(NKcell_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_NKcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
NKcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_NKcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
NKcell_simp <- rbind(NKcell_D2v0_simp, NKcell_D8v0_simp, NKcell_D8v2_simp)
write.table(NKcell_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_NKcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
Bcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_all <- rbind(Bcell_D2v0, Bcell_D8v0, Bcell_D8v2)
write.table(Bcell_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Bcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
Bcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_Bcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
Bcell_simp <- rbind(Bcell_D2v0_simp, Bcell_D8v0_simp, Bcell_D8v2_simp)
write.table(Bcell_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_Bcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
gdTcell_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_all <- rbind(gdTcell_D2v0, gdTcell_D8v0, gdTcell_D8v2)
write.table(gdTcell_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD2neg_GD_Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
gdTcell_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
gdTcell_simp <- rbind(gdTcell_D2v0_simp, gdTcell_D8v0_simp, gdTcell_D8v2_simp)
write.table(gdTcell_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD2neg_GD_Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
cDCs_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
cDCs_all <- cDCs_D2v0
write.table(cDCs_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_cDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
cDCs_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_cDCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
cDCs_simp <- cDCs_D2v0_simp
write.table(cDCs_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_cDCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
pDCs_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_all <- rbind(pDCs_D2v0, pDCs_D8v0, pDCs_D8v2)
write.table(pDCs_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_pDCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
pDCs_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_pDCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_pDCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_pDCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
pDCs_simp <- rbind(pDCs_D2v0_simp, pDCs_D8v0_simp, pDCs_D8v2_simp)
write.table(pDCs_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_pDCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
CD4Tcells_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_all <- rbind(CD4Tcells_D2v0, CD4Tcells_D8v0, CD4Tcells_D8v2)
write.table(CD4Tcells_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD4Tcells_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
CD4Tcells_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
CD4Tcells_simp <- rbind(CD4Tcells_D2v0_simp, CD4Tcells_D8v0_simp, CD4Tcells_D8v2_simp)
write.table(CD4Tcells_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_CD4Tcells_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)
ASCs_D2v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_D8v0 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_D8v2 <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_all <- rbind(ASCs_D2v0, ASCs_D8v0, ASCs_D8v2)
write.table(ASCs_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_ASCs_pos_enrichGO_2026_04_14.txt", sep = "\t", row.names = FALSE)
ASCs_D2v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D2v0_ASCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_D8v0_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v0_ASCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_D8v2_simp <- read.table("/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/D8v2_ASCs_pos_enrichGO_simplified_2026_04_14.txt", header=TRUE, sep="\t")
ASCs_simp <- rbind(ASCs_D2v0_simp, ASCs_D8v0_simp, ASCs_D8v2_simp)
write.table(ASCs_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_ASCs_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)

rbind_all <- rbind(Monocyte_all, NKcell_all, Bcell_all, gdTcell_all, cDCs_all, pDCs_all, CD4Tcells_all, ASCs_all)
write.table(rbind_all, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_celltypes_pos_enrichGO_all_2026_04_14.txt", sep = "\t", row.names = FALSE)
rbind_simp <- rbind(Monocyte_simp, NKcell_simp, Bcell_simp, gdTcell_simp, cDCs_simp, pDCs_simp, CD4Tcells_simp, ASCs_simp)
write.table(rbind_simp, file="/scRNAseq/Sal_5pigs_2026/GO/between_timepoints/All_celltypes_pos_enrichGO_simplified_2026_04_14.txt", sep = "\t", row.names = FALSE)

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