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
pigrnaatlas<-read.table("/Annotation_files/pig_rna_atlas.tsv", sep = "\t", header = TRUE)
#keep columns "Gene", "Gene.description"
pigrnaatlas2<-pigrnaatlas[,c("Gene","Ensembl","Gene.description")]
#read in the orthologs
ORG<- read.csv(file ="/Annotation_files/PigToHuman_GeneOrthos_v11_1_97_scGenes.csv", header = T,row.names=1)
colnames(ORG)[11] <- "Gene"
# Replace "_" with "-" in the Gene column
ORG$Gene <- gsub("_", "-", ORG$Gene)
# For duplicates in the Gene column, keep the row with the highest X.id..query.gene.identical.to.target.Human.gene
#ORG2 <- ORG %>% group_by(Gene) %>% filter(X.id..query.gene.identical.to.target.Human.gene == max(X.id..query.gene.identical.to.target.Human.gene)) %>% ungroup()
#subset the columns needed for annotation
ORG2.subset <- ORG %>% dplyr::select(c(1,11,7,8,10))
ORG2.subset$Gene <- as.character(ORG2.subset$Gene)

Monocytes_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Monocytes_D2v0_2<-left_join(Monocytes_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(Monocytes_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

Monocytes_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Monocytes_D8v2_2<-left_join(Monocytes_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(Monocytes_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

Monocytes_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Monocytes_D8v0_2<-left_join(Monocytes_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(Monocytes_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
#combine
Monocytes_day <- rbind(Monocytes_D2v0_2, Monocytes_D8v2_2, Monocytes_D8v0_2)
Monocytes_day<-write.table(Monocytes_day, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

NKcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
NKcells_D2v0_2<-left_join(NKcells_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(NKcells_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

NKcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
NKcells_D8v2_2<-left_join(NKcells_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(NKcells_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

NKcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
NKcells_D8v0_2<-left_join(NKcells_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(NKcells_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
#combine
NKcells_day <- rbind(NKcells_D2v0_2, NKcells_D8v2_2, NKcells_D8v0_2)
NKcells_day<-write.table(NKcells_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

Bcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Bcells_D2v0_2<-left_join(Bcells_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(Bcells_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

Bcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Bcells_D8v2_2<-left_join(Bcells_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(Bcells_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

Bcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
Bcells_D8v0_2<-left_join(Bcells_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(Bcells_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
#combine
Bcells_day <- rbind(Bcells_D2v0_2, Bcells_D8v2_2, Bcells_D8v0_2)
Bcells_day <- write.table(Bcells_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

CD2neg_GD_Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD2neg_GD_Tcells_D2v0_2<-left_join(CD2neg_GD_Tcells_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(CD2neg_GD_Tcells_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

CD2neg_GD_Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD2neg_GD_Tcells_D8v2_2<-left_join(CD2neg_GD_Tcells_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(CD2neg_GD_Tcells_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

CD2neg_GD_Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD2neg_GD_Tcells_D8v0_2<-left_join(CD2neg_GD_Tcells_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(CD2neg_GD_Tcells_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
#combine
CD2neg_GD_Tcells_day <- rbind(CD2neg_GD_Tcells_D2v0_2, CD2neg_GD_Tcells_D8v2_2, CD2neg_GD_Tcells_D8v0_2)
CD2neg_GD_Tcells_day <- write.table(CD2neg_GD_Tcells_day,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

ASCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
ASCs_D2v0_2<-left_join(ASCs_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(ASCs_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

ASCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
ASCs_D8v2_2<-left_join(ASCs_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(ASCs_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

ASCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
ASCs_D8v0_2<-left_join(ASCs_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(ASCs_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
#combine
ASCs_day <- rbind(ASCs_D2v0_2, ASCs_D8v2_2, ASCs_D8v0_2)
ASCs_day <- write.table(ASCs_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

pDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
pDCs_D2v0_2<-left_join(pDCs_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(pDCs_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

pDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
pDCs_D8v2_2<-left_join(pDCs_D8v2,ORG2.subset,by=c("gene"="Gene"))
if(nrow(pDCs_D8v2_2) > 0) {
write.table(pDCs_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(pDCs_D8v2,pDCs_D8v2_2)
}

pDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
pDCs_D8v0_2<-left_join(pDCs_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(pDCs_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
#combine
# Check if the data frames exist before combining
if (exists("pDCs_D2v0_2") && exists("pDCs_D8v2_2") && exists("pDCs_D8v0_2")) {
  pDCs_day <- rbind(pDCs_D2v0_2, pDCs_D8v2_2, pDCs_D8v0_2)
  write.table(pDCs_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t", row.names=FALSE)
} else {
  warning("One or more data frames are missing. Skipping rbind.")
}

cDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
if (exists("cDCs_D8v2")) {
cDCs_D8v2_2<-left_join(cDCs_D8v2,ORG2.subset,by=c("gene"="Gene"))
}
if (exists("cDCs_D8v2_2 ") && nrow(cDCs_D8v2_2) > 0) {
write.table(cDCs_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
}

cDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
if (exists("cDCs_D8v0")) {
cDCs_D8v0$"gene" <- as.character(cDCs_D8v0$"gene")
cDCs_D8v0_2<-left_join(cDCs_D8v0,ORG2.subset,by=c("gene"="Gene"))
}
if (exists("cDCs_D8v0_2") && nrow(cDCs_D8v0_2) > 0) {
write.table(cDCs_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
}
#combine
if (exists("cDCs_D8v2_2") && exists("cDCs_D8v0_2")) {
cDCs_day <- rbind(cDCs_D8v2_2, cDCs_D8v0_2)
cDCs_day <- write.table(cDCs_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
} else {
  warning("One or more cDCs data frames are missing. Skipping rbind.")
}

CD4Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD4Tcells_D2v0_2<-left_join(CD4Tcells_D2v0,ORG2.subset,by=c("gene"="Gene"))
write.table(CD4Tcells_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

CD4Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD4Tcells_D8v2_2<-left_join(CD4Tcells_D8v2,ORG2.subset,by=c("gene"="Gene"))
write.table(CD4Tcells_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

CD4Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
CD4Tcells_D8v0_2<-left_join(CD4Tcells_D8v0,ORG2.subset,by=c("gene"="Gene"))
write.table(CD4Tcells_D8v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
#combine
CD4Tcells_day <- rbind(CD4Tcells_D2v0_2, CD4Tcells_D8v2_2, CD4Tcells_D8v0_2)
CD4Tcells_day <- write.table(CD4Tcells_day, file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_allComparisons_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)


######################### Merge all significant DE gene files and add cell type and comparison columns #########################
Monocytes_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

Monocytes_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

Monocytes_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

NKcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

NKcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

NKcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

Bcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

Bcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

Bcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D2v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D8v2 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D8v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

ASCs_D2v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

ASCs_D8v2 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

ASCs_D8v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

pDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

pDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

pDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

cDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)


cDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

cDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD4Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD4Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

CD4Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

# Combine all data frames into one
all_sig <- rbind(ASCs_D2v0, ASCs_D8v0, ASCs_D8v2,Bcells_D2v0, Bcells_D8v0, Bcells_D8v2,CD2neg_GD_Tcells_D2v0, CD2neg_GD_Tcells_D8v0, CD2neg_GD_Tcells_D8v2,CD4Tcells_D2v0, CD4Tcells_D8v0, CD4Tcells_D8v2,cDCs_D2v0, cDCs_D8v0, cDCs_D8v2,Monocytes_D2v0, Monocytes_D8v0, Monocytes_D8v2,NKcells_D2v0, NKcells_D8v0, NKcells_D8v2,pDCs_D2v0, pDCs_D8v0, pDCs_D8v2)
# Write the summary to a file
write.table(all_sig, "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Sal_FindMarkers_pairwise_SIG_DEGs_allCelltypes.txt", sep = "\t", quote = FALSE, row.names = TRUE)

######################### Get number of DEGs per cell type and comparison by Log2FC #########################
Monocytes_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Monocytes D2v0
Monocytes_D2v0_posLogFC <- Monocytes_D2v0[which(Monocytes_D2v0$avg_log2FC > 0),]
Monocytes_D2v0_posLogFC_count <- nrow(Monocytes_D2v0_posLogFC)
#get the number of genes with negative logFC in Monocytes D2v0
Monocytes_D2v0_negLogFC <- Monocytes_D2v0[which(Monocytes_D2v0$avg_log2FC < 0),]
Monocytes_D2v0_negLogFC_count <- nrow(Monocytes_D2v0_negLogFC)

Monocytes_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Monocytes D8v2
Monocytes_D8v2_posLogFC <- Monocytes_D8v2[which(Monocytes_D8v2$avg_log2FC > 0),]
Monocytes_D8v2_posLogFC_count <- nrow(Monocytes_D8v2_posLogFC)
#get the number of genes with negative logFC in Monocytes D8v2
Monocytes_D8v2_negLogFC <- Monocytes_D8v2[which(Monocytes_D8v2$avg_log2FC < 0),]
Monocytes_D8v2_negLogFC_count <- nrow(Monocytes_D8v2_negLogFC)

Monocytes_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Monocytes D8v0
Monocytes_D8v0_posLogFC <- Monocytes_D8v0[which(Monocytes_D8v0$avg_log2FC > 0),]
Monocytes_D8v0_posLogFC_count <- nrow(Monocytes_D8v0_posLogFC)
#get the number of genes with negative logFC in Monocytes D8v0
Monocytes_D8v0_negLogFC <- Monocytes_D8v0[which(Monocytes_D8v0$avg_log2FC < 0),]
Monocytes_D8v0_negLogFC_count <- nrow(Monocytes_D8v0_negLogFC)

NKcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in NKcells D2v0
NKcells_D2v0_posLogFC <- NKcells_D2v0[which(NKcells_D2v0$avg_log2FC > 0),]
NKcells_D2v0_posLogFC_count <- nrow(NKcells_D2v0_posLogFC)
#get the number of genes with negative logFC in NKcells D2v0
NKcells_D2v0_negLogFC <- NKcells_D2v0[which(NKcells_D2v0$avg_log2FC < 0),]
NKcells_D2v0_negLogFC_count <- nrow(NKcells_D2v0_negLogFC)

NKcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in NKcells D8v2
NKcells_D8v2_posLogFC <- NKcells_D8v2[which(NKcells_D8v2$avg_log2FC > 0),]
NKcells_D8v2_posLogFC_count <- nrow(NKcells_D8v2_posLogFC)
#get the number of genes with negative logFC in NKcells D8v2
NKcells_D8v2_negLogFC <- NKcells_D8v2[which(NKcells_D8v2$avg_log2FC < 0),]
NKcells_D8v2_negLogFC_count <- nrow(NKcells_D8v2_negLogFC)

NKcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in NKcells D8v0
NKcells_D8v0_posLogFC <- NKcells_D8v0[which(NKcells_D8v0$avg_log2FC > 0),]
NKcells_D8v0_posLogFC_count <- nrow(NKcells_D8v0_posLogFC)
#get the number of genes with negative logFC in NKcells D8v0
NKcells_D8v0_negLogFC <- NKcells_D8v0[which(NKcells_D8v0$avg_log2FC < 0),]
NKcells_D8v0_negLogFC_count <- nrow(NKcells_D8v0_negLogFC)

Bcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Bcells D2v0
Bcells_D2v0_posLogFC <- Bcells_D2v0[which(Bcells_D2v0$avg_log2FC > 0),]
Bcells_D2v0_posLogFC_count <- nrow(Bcells_D2v0_posLogFC)
#get the number of genes with negative logFC in Bcells D2v0
Bcells_D2v0_negLogFC <- Bcells_D2v0[which(Bcells_D2v0$avg_log2FC < 0),]
Bcells_D2v0_negLogFC_count <- nrow(Bcells_D2v0_negLogFC)

Bcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Bcells D8v2
Bcells_D8v2_posLogFC <- Bcells_D8v2[which(Bcells_D8v2$avg_log2FC > 0),]
Bcells_D8v2_posLogFC_count <- nrow(Bcells_D8v2_posLogFC)
#get the number of genes with negative logFC in Bcells D8v2
Bcells_D8v2_negLogFC <- Bcells_D8v2[which(Bcells_D8v2$avg_log2FC < 0),]
Bcells_D8v2_negLogFC_count <- nrow(Bcells_D8v2_negLogFC)

Bcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in Bcells D8v0
Bcells_D8v0_posLogFC <- Bcells_D8v0[which(Bcells_D8v0$avg_log2FC > 0),]
Bcells_D8v0_posLogFC_count <- nrow(Bcells_D8v0_posLogFC)
#get the number of genes with negative logFC in Bcells D8v0
Bcells_D8v0_negLogFC <- Bcells_D8v0[which(Bcells_D8v0$avg_log2FC < 0),]
Bcells_D8v0_negLogFC_count <- nrow(Bcells_D8v0_negLogFC)

CD2neg_GD_Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in GD Tcells D2v0
CD2neg_GD_Tcells_D2v0_posLogFC <- CD2neg_GD_Tcells_D2v0[which(CD2neg_GD_Tcells_D2v0$avg_log2FC > 0),]
CD2neg_GD_Tcells_D2v0_posLogFC_count <- nrow(CD2neg_GD_Tcells_D2v0_posLogFC)
#get the number of genes with negative logFC in GD Tcells D2v0
CD2neg_GD_Tcells_D2v0_negLogFC <- CD2neg_GD_Tcells_D2v0[which(CD2neg_GD_Tcells_D2v0$avg_log2FC < 0),]
CD2neg_GD_Tcells_D2v0_negLogFC_count <- nrow(CD2neg_GD_Tcells_D2v0_negLogFC)

CD2neg_GD_Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in GD Tcells D8v2
CD2neg_GD_Tcells_D8v2_posLogFC <- CD2neg_GD_Tcells_D8v2[which(CD2neg_GD_Tcells_D8v2$avg_log2FC > 0),]
CD2neg_GD_Tcells_D8v2_posLogFC_count <- nrow(CD2neg_GD_Tcells_D8v2_posLogFC)
#get the number of genes with negative logFC in GD Tcells D8v2
CD2neg_GD_Tcells_D8v2_negLogFC <- CD2neg_GD_Tcells_D8v2[which(CD2neg_GD_Tcells_D8v2$avg_log2FC < 0),]
CD2neg_GD_Tcells_D8v2_negLogFC_count <- nrow(CD2neg_GD_Tcells_D8v2_negLogFC)

CD2neg_GD_Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in GD Tcells D8v0
CD2neg_GD_Tcells_D8v0_posLogFC <- CD2neg_GD_Tcells_D8v0[which(CD2neg_GD_Tcells_D8v0$avg_log2FC > 0),]
CD2neg_GD_Tcells_D8v0_posLogFC_count <- nrow(CD2neg_GD_Tcells_D8v0_posLogFC)
#get the number of genes with negative logFC in GD Tcells D8v0
CD2neg_GD_Tcells_D8v0_negLogFC <- CD2neg_GD_Tcells_D8v0[which(CD2neg_GD_Tcells_D8v0$avg_log2FC < 0),]
CD2neg_GD_Tcells_D8v0_negLogFC_count <- nrow(CD2neg_GD_Tcells_D8v0_negLogFC)

ASCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in ASCs D2v0
ASCs_D2v0_posLogFC <- ASCs_D2v0[which(ASCs_D2v0$avg_log2FC > 0),]
ASCs_D2v0_posLogFC_count <- nrow(ASCs_D2v0_posLogFC)
#get the number of genes with negative logFC in ASCs D2v0
ASCs_D2v0_negLogFC <- ASCs_D2v0[which(ASCs_D2v0$avg_log2FC < 0),]
ASCs_D2v0_negLogFC_count <- nrow(ASCs_D2v0_negLogFC)

ASCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in ASCs D8v2
ASCs_D8v2_posLogFC <- ASCs_D8v2[which(ASCs_D8v2$avg_log2FC > 0),]
ASCs_D8v2_posLogFC_count <- nrow(ASCs_D8v2_posLogFC)
#get the number of genes with negative logFC in ASCs D8v2
ASCs_D8v2_negLogFC <- ASCs_D8v2[which(ASCs_D8v2$avg_log2FC < 0),]
ASCs_D8v2_negLogFC_count <- nrow(ASCs_D8v2_negLogFC)

ASCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in ASCs D8v0
ASCs_D8v0_posLogFC <- ASCs_D8v0[which(ASCs_D8v0$avg_log2FC > 0),]
ASCs_D8v0_posLogFC_count <- nrow(ASCs_D8v0_posLogFC)
#get the number of genes with negative logFC in ASCs D8v0
ASCs_D8v0_negLogFC <- ASCs_D8v0[which(ASCs_D8v0$avg_log2FC < 0),]
ASCs_D8v0_negLogFC_count <- nrow(ASCs_D8v0_negLogFC)

pDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in pDCs D2v0
pDCs_D2v0_posLogFC <- pDCs_D2v0[which(pDCs_D2v0$avg_log2FC > 0),]
pDCs_D2v0_posLogFC_count <- nrow(pDCs_D2v0_posLogFC)
#get the number of genes with negative logFC in pDCs D2v0
pDCs_D2v0_negLogFC <- pDCs_D2v0[which(pDCs_D2v0$avg_log2FC < 0),]
pDCs_D2v0_negLogFC_count <- nrow(pDCs_D2v0_negLogFC)

pDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in pDCs D8v2
pDCs_D8v2_posLogFC <- pDCs_D8v2[which(pDCs_D8v2$avg_log2FC > 0),]
pDCs_D8v2_posLogFC_count <- nrow(pDCs_D8v2_posLogFC)
#get the number of genes with negative logFC in pDCs D8v2
pDCs_D8v2_negLogFC <- pDCs_D8v2[which(pDCs_D8v2$avg_log2FC < 0),]
pDCs_D8v2_negLogFC_count <- nrow(pDCs_D8v2_negLogFC)

pDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in pDCs D8v0
pDCs_D8v0_posLogFC <- pDCs_D8v0[which(pDCs_D8v0$avg_log2FC > 0),]
pDCs_D8v0_posLogFC_count <- nrow(pDCs_D8v0_posLogFC)
#get the number of genes with negative logFC in pDCs D8v0
pDCs_D8v0_negLogFC <- pDCs_D8v0[which(pDCs_D8v0$avg_log2FC < 0),]
pDCs_D8v0_negLogFC_count <- nrow(pDCs_D8v0_negLogFC)

cDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in cDCs D2v0
cDCs_D2v0_posLogFC <- cDCs_D2v0[which(cDCs_D2v0$avg_log2FC > 0),]
cDCs_D2v0_posLogFC_count <- nrow(cDCs_D2v0_posLogFC)
#get the number of genes with negative logFC in cDCs D2v0
cDCs_D2v0_negLogFC <- cDCs_D2v0[which(cDCs_D2v0$avg_log2FC < 0),]
cDCs_D2v0_negLogFC_count <- nrow(cDCs_D2v0_negLogFC)

cDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in cDCs D8v2
cDCs_D8v2_posLogFC <- cDCs_D8v2[which(cDCs_D8v2$avg_log2FC > 0),]
cDCs_D8v2_posLogFC_count <- nrow(cDCs_D8v2_posLogFC)
#get the number of genes with negative logFC in cDCs D8v2
cDCs_D8v2_negLogFC <- cDCs_D8v2[which(cDCs_D8v2$avg_log2FC < 0),]
cDCs_D8v2_negLogFC_count <- nrow(cDCs_D8v2_negLogFC)

cDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in cDCs D8v0
cDCs_D8v0_posLogFC <- cDCs_D8v0[which(cDCs_D8v0$avg_log2FC > 0),]
cDCs_D8v0_posLogFC_count <- nrow(cDCs_D8v0_posLogFC)
#get the number of genes with negative logFC in cDCs D8v0
cDCs_D8v0_negLogFC <- cDCs_D8v0[which(cDCs_D8v0$avg_log2FC < 0),]
cDCs_D8v0_negLogFC_count <- nrow(cDCs_D8v0_negLogFC)

CD4Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in CD4Tcells D2v0
CD4Tcells_D2v0_posLogFC <- CD4Tcells_D2v0[which(CD4Tcells_D2v0$avg_log2FC > 0),]
CD4Tcells_D2v0_posLogFC_count <- nrow(CD4Tcells_D2v0_posLogFC)
#get the number of genes with negative logFC in CD4Tcells D2v0
CD4Tcells_D2v0_negLogFC <- CD4Tcells_D2v0[which(CD4Tcells_D2v0$avg_log2FC < 0),]
CD4Tcells_D2v0_negLogFC_count <- nrow(CD4Tcells_D2v0_negLogFC)

CD4Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in CD4Tcells D8v2
CD4Tcells_D8v2_posLogFC <- CD4Tcells_D8v2[which(CD4Tcells_D8v2$avg_log2FC > 0),]
CD4Tcells_D8v2_posLogFC_count <- nrow(CD4Tcells_D8v2_posLogFC)
#get the number of genes with negative logFC in CD4Tcells D8v2
CD4Tcells_D8v2_negLogFC <- CD4Tcells_D8v2[which(CD4Tcells_D8v2$avg_log2FC < 0),]
CD4Tcells_D8v2_negLogFC_count <- nrow(CD4Tcells_D8v2_negLogFC)

CD4Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
#get the number of genes with positive logFC in CD4Tcells D8v0
CD4Tcells_D8v0_posLogFC <- CD4Tcells_D8v0[which(CD4Tcells_D8v0$avg_log2FC > 0),]
CD4Tcells_D8v0_posLogFC_count <- nrow(CD4Tcells_D8v0_posLogFC)
#get the number of genes with negative logFC in CD4Tcells D8v0
CD4Tcells_D8v0_negLogFC <- CD4Tcells_D8v0[which(CD4Tcells_D8v0$avg_log2FC < 0),]
CD4Tcells_D8v0_negLogFC_count <- nrow(CD4Tcells_D8v0_negLogFC)

# Combine all _posLogFC & _negLogFC vectors 
All_GO_DEGs <- rbind(Monocytes_D2v0_posLogFC_count,Monocytes_D2v0_negLogFC_count, Monocytes_D8v2_posLogFC_count,Monocytes_D8v2_negLogFC_count,Monocytes_D8v0_posLogFC_count,Monocytes_D8v0_negLogFC_count,NKcells_D2v0_posLogFC_count, NKcells_D2v0_negLogFC_count, NKcells_D8v2_posLogFC_count,  NKcells_D8v2_negLogFC_count, NKcells_D8v0_posLogFC_count,NKcells_D8v0_negLogFC_count,Bcells_D2v0_posLogFC_count,Bcells_D2v0_negLogFC_count, Bcells_D8v2_posLogFC_count, Bcells_D8v2_negLogFC_count, Bcells_D8v0_posLogFC_count, Bcells_D8v0_negLogFC_count,CD2neg_GD_Tcells_D2v0_posLogFC_count,CD2neg_GD_Tcells_D2v0_negLogFC_count,CD2neg_GD_Tcells_D8v2_posLogFC_count,CD2neg_GD_Tcells_D8v2_negLogFC_count,CD2neg_GD_Tcells_D8v0_posLogFC_count,  CD2neg_GD_Tcells_D8v0_negLogFC_count, ASCs_D2v0_posLogFC_count,ASCs_D2v0_negLogFC_count, ASCs_D8v2_posLogFC_count, ASCs_D8v2_negLogFC_count,ASCs_D8v0_posLogFC_count, ASCs_D8v0_negLogFC_count,pDCs_D2v0_posLogFC_count, pDCs_D8v2_posLogFC_count, pDCs_D2v0_posLogFC_count,pDCs_D2v0_negLogFC_count, pDCs_D8v2_posLogFC_count,pDCs_D8v2_negLogFC_count,pDCs_D8v0_posLogFC_count,pDCs_D8v0_negLogFC_count,cDCs_D2v0_posLogFC_count, cDCs_D2v0_negLogFC_count,cDCs_D8v2_posLogFC_count,cDCs_D8v2_negLogFC_count, cDCs_D8v0_posLogFC_count,cDCs_D8v0_negLogFC_count,CD4Tcells_D2v0_posLogFC_count,CD4Tcells_D2v0_negLogFC_count,CD4Tcells_D8v2_posLogFC_count,CD4Tcells_D8v2_negLogFC_count,CD4Tcells_D8v0_posLogFC_count,CD4Tcells_D8v0_negLogFC_count)

All_GO_DEGs
#save as a text file keep row names

write.table(All_GO_DEGs, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Number_of_DEGs_PosLogFC_NegLogFC.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)
