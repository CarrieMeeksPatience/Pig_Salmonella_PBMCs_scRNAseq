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

sal<-readRDS("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")

DefaultAssay(sal)<-"RNA"
#read in the pig rna atlas
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

Idents(sal)<-"celltypes_subclustres02"
print("table(sal$celltypes_subclustres02)")
table(sal$celltypes_subclustres02)

Monocytes<-subset(sal, idents = c("Monocytes"))
#Get names of all genes in Monocytes 
Monocytes_bkg<-Monocytes[["RNA"]]$counts
#remove non expressed genes
Monocytes_bkg<-Monocytes_bkg[rowSums(Monocytes_bkg) > 0, ]
Monocytes_bkg<-as.data.frame(Monocytes_bkg)
Monocytes_bkg<-rownames(Monocytes)
Monocytes_bkg<-as.data.frame(Monocytes_bkg)
colnames(Monocytes_bkg)<-"Gene"
#Add human geneID to the Monocytes_bkg
Monocytes_bkg<-left_join(Monocytes_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
Monocytes_bkg<-Monocytes_bkg[!duplicated(Monocytes_bkg$Gene), ]
write.table(Monocytes_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(Monocytes)<-"time_point"
print("table(Idents(Monocytes))")
table(Idents(Monocytes))
#Run FindMarkers on Monocytes
Monocytes_D2v0<-FindMarkers(object = Monocytes, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Monocytes_D2v0$"pct.1-pct.2" <- Monocytes_D2v0$pct.1 - Monocytes_D2v0$pct.2
Monocytes_D2v0$"gene" <- rownames(Monocytes_D2v0)
Monocytes_D2v0$Celltype <- "Monocytes"
Monocytes_D2v0$Comparison <- "2DPI vs 0DPI"

write.csv(Monocytes_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers.csv")
#save as txt file
write.table(Monocytes_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
Monocytes_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Monocytes_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Monocytes_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")
Monocytes_D8v2<-FindMarkers(Monocytes, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Monocytes_D8v2$"pct.1-pct.2" <- Monocytes_D8v2$pct.1 - Monocytes_D8v2$pct.2
Monocytes_D8v2$"gene" <- rownames(Monocytes_D8v2)
Monocytes_D8v2$Celltype <- "Monocytes"
Monocytes_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(Monocytes_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers.csv")
#save as txt file
write.table(Monocytes_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(Monocytes_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Monocytes_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
Monocytes_D8v0<-FindMarkers(Monocytes, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Monocytes_D8v0$"pct.1-pct.2" <- Monocytes_D8v0$pct.1 - Monocytes_D8v0$pct.2
Monocytes_D8v0$"gene" <- rownames(Monocytes_D8v0)
Monocytes_D8v0$Celltype <- "Monocytes"
Monocytes_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(Monocytes_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers.csv")
#save as txt file
write.table(Monocytes_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
Monocytes_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Monocytes_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Monocytes_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)

##NK cells
print("Running DE on NK cells")
NKcells<-subset(sal, idents = c("NK cells"))
#Get names of all genes in NKcells 
NKcells_bkg<-NKcells[["RNA"]]$counts
#remove non expressed genes
NKcells_bkg<-NKcells_bkg[rowSums(NKcells_bkg) > 0, ]
NKcells_bkg<-as.data.frame(NKcells_bkg)
NKcells_bkg<-rownames(NKcells)
NKcells_bkg<-as.data.frame(NKcells_bkg)
colnames(NKcells_bkg)<-"Gene"
#Add human geneID to the NKcells_bkg
NKcells_bkg<-left_join(NKcells_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
NKcells_bkg<-NKcells_bkg[!duplicated(NKcells_bkg$Gene), ]
write.table(NKcells_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(NKcells)<-"time_point"
print("table(Idents(NKcells))")
table(Idents(NKcells))
#Run FindMarkers on NKcells
NKcells_D2v0<-FindMarkers(object = NKcells, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
NKcells_D2v0$"pct.1-pct.2" <- NKcells_D2v0$pct.1 - NKcells_D2v0$pct.2
NKcells_D2v0$"gene" <- rownames(NKcells_D2v0)
NKcells_D2v0$Celltype <- "NK cells"
NKcells_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(NKcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers.csv")
#save as txt file
write.table(NKcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
NKcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(NKcells_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(NKcells_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")
NKcells_D8v2<-FindMarkers(NKcells, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
NKcells_D8v2$"pct.1-pct.2" <- NKcells_D8v2$pct.1 - NKcells_D8v2$pct.2
NKcells_D8v2$"gene" <- rownames(NKcells_D8v2)
NKcells_D8v2$Celltype <- "NK cells"
NKcells_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(NKcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers.csv")
#save as txt file
write.table(NKcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(NKcells_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(NKcells_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
NKcells_D8v0<-FindMarkers(NKcells, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
NKcells_D8v0$"pct.1-pct.2" <- NKcells_D8v0$pct.1 - NKcells_D8v0$pct.2
NKcells_D8v0$"gene" <- rownames(NKcells_D8v0)
NKcells_D8v0$Celltype <- "NK cells"
NKcells_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(NKcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers.csv")
#save as txt file
write.table(NKcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
NKcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(NKcells_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(NKcells_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)


## B cells
print("running DE analysis on B cells")
Bcells<-subset(sal, idents = c("B cells"))
#Get names of all genes in Bcells 
Bcells_bkg<-Bcells[["RNA"]]$counts
#remove non expressed genes
Bcells_bkg<-Bcells_bkg[rowSums(Bcells_bkg) > 0, ]
Bcells_bkg<-as.data.frame(Bcells_bkg)
Bcells_bkg<-rownames(Bcells)
Bcells_bkg<-as.data.frame(Bcells_bkg)
colnames(Bcells_bkg)<-"Gene"
#Add human geneID to the Bcells_bkg
Bcells_bkg<-left_join(Bcells_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
Bcells_bkg<-Bcells_bkg[!duplicated(Bcells_bkg$Gene), ]
write.table(Bcells_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(Bcells)<-"time_point"
print("table(Idents(Bcells))")
table(Idents(Bcells))
#Run FindMarkers on Bcells
Bcells_D2v0<-FindMarkers(object = Bcells, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Bcells_D2v0$"pct.1-pct.2" <- Bcells_D2v0$pct.1 - Bcells_D2v0$pct.2
Bcells_D2v0$"gene" <- rownames(Bcells_D2v0)
Bcells_D2v0$Celltype <- "B cells"
Bcells_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(Bcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers.csv")
#save as txt file
write.table(Bcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
Bcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Bcells_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Bcells_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")
Bcells_D8v2<-FindMarkers(Bcells, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Bcells_D8v2$"pct.1-pct.2" <- Bcells_D8v2$pct.1 - Bcells_D8v2$pct.2
Bcells_D8v2$"gene" <- rownames(Bcells_D8v2)
Bcells_D8v2$Celltype <- "B cells"
Bcells_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(Bcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers.csv")
#save as txt file
write.table(Bcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(Bcells_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Bcells_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
Bcells_D8v0<-FindMarkers(Bcells, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Bcells_D8v0$"pct.1-pct.2" <- Bcells_D8v0$pct.1 - Bcells_D8v0$pct.2
Bcells_D8v0$"gene" <- rownames(Bcells_D8v0)
Bcells_D8v0$Celltype <- "B cells"
Bcells_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(Bcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers.csv")
#save as txt file
write.table(Bcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
Bcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Bcells_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Bcells_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)




## GD T-cells
print("running analysis on  GD T-cells")
CD2neg_GD_Tcells<-subset(sal, idents = c("CD2- gd T cells"))
#Get names of all genes in CD2neg_GD_Tcells 
CD2neg_GD_Tcells_bkg<-CD2neg_GD_Tcells[["RNA"]]$counts
#remove non expressed genes
CD2neg_GD_Tcells_bkg<-CD2neg_GD_Tcells_bkg[rowSums(CD2neg_GD_Tcells_bkg) > 0, ]
CD2neg_GD_Tcells_bkg<-as.data.frame(CD2neg_GD_Tcells_bkg)
CD2neg_GD_Tcells_bkg<-rownames(CD2neg_GD_Tcells)
CD2neg_GD_Tcells_bkg<-as.data.frame(CD2neg_GD_Tcells_bkg)
colnames(CD2neg_GD_Tcells_bkg)<-"Gene"
#Add human geneID to the CD2neg_GD_Tcells_bkg
CD2neg_GD_Tcells_bkg<-left_join(CD2neg_GD_Tcells_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
CD2neg_GD_Tcells_bkg<-CD2neg_GD_Tcells_bkg[!duplicated(CD2neg_GD_Tcells_bkg$Gene), ]
write.table(CD2neg_GD_Tcells_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(CD2neg_GD_Tcells)<-"time_point"
print("table(Idents(CD2neg_GD_Tcells))")
table(Idents(CD2neg_GD_Tcells))
#Run FindMarkers on CD2neg_GD_Tcells
CD2neg_GD_Tcells_D2v0<-FindMarkers(object = CD2neg_GD_Tcells, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD2neg_GD_Tcells_D2v0$"pct.1-pct.2" <- CD2neg_GD_Tcells_D2v0$pct.1 - CD2neg_GD_Tcells_D2v0$pct.2
CD2neg_GD_Tcells_D2v0$"gene" <- rownames(CD2neg_GD_Tcells_D2v0)
CD2neg_GD_Tcells_D2v0$Celltype <- "CD2- gd T cells"
CD2neg_GD_Tcells_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(CD2neg_GD_Tcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers.csv")
#save as txt file
write.table(CD2neg_GD_Tcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
CD2neg_GD_Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(CD2neg_GD_Tcells_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD2neg_GD_Tcells_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")

CD2neg_GD_Tcells_D8v2<-FindMarkers(CD2neg_GD_Tcells, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD2neg_GD_Tcells_D8v2$"pct.1-pct.2" <- CD2neg_GD_Tcells_D8v2$pct.1 - CD2neg_GD_Tcells_D8v2$pct.2
CD2neg_GD_Tcells_D8v2$"gene" <- rownames(CD2neg_GD_Tcells_D8v2)
CD2neg_GD_Tcells_D8v2$Celltype <- "CD2- gd T cells"
CD2neg_GD_Tcells_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(CD2neg_GD_Tcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers.csv")
#save as txt file
write.table(CD2neg_GD_Tcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(CD2neg_GD_Tcells_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD2neg_GD_Tcells_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
CD2neg_GD_Tcells_D8v0<-FindMarkers(CD2neg_GD_Tcells, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD2neg_GD_Tcells_D8v0$"pct.1-pct.2" <- CD2neg_GD_Tcells_D8v0$pct.1 - CD2neg_GD_Tcells_D8v0$pct.2
CD2neg_GD_Tcells_D8v0$"gene" <- rownames(CD2neg_GD_Tcells_D8v0)
CD2neg_GD_Tcells_D8v0$Celltype <- "CD2- gd T cells"
CD2neg_GD_Tcells_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(CD2neg_GD_Tcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers.csv")
#save as txt file
write.table(CD2neg_GD_Tcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
CD2neg_GD_Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(CD2neg_GD_Tcells_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD2neg_GD_Tcells_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)



## ASCs
print("running analysis on ASCs")
ASCs<-subset(sal, idents = c("ASCs"))
#Get names of all genes in ASCs 
ASCs_bkg<-ASCs[["RNA"]]$counts
#remove non expressed genes
ASCs_bkg<-ASCs_bkg[rowSums(ASCs_bkg) > 0, ]
ASCs_bkg<-as.data.frame(ASCs_bkg)
ASCs_bkg<-rownames(ASCs)
ASCs_bkg<-as.data.frame(ASCs_bkg)
colnames(ASCs_bkg)<-"Gene"
#Add human geneID to the ASCs_bkg
ASCs_bkg<-left_join(ASCs_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
ASCs_bkg<-ASCs_bkg[!duplicated(ASCs_bkg$Gene), ]
write.table(ASCs_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(ASCs)<-"time_point"
print("table(Idents(ASCs))")
table(Idents(ASCs))
#Run FindMarkers on ASCs
ASCs_D2v0<-FindMarkers(object = ASCs, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
ASCs_D2v0$"pct.1-pct.2" <- ASCs_D2v0$pct.1 - ASCs_D2v0$pct.2
ASCs_D2v0$"gene" <- rownames(ASCs_D2v0)
ASCs_D2v0$Celltype <- "ASCs"
ASCs_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(ASCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers.csv")
#save as txt file
write.table(ASCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
ASCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(ASCs_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(ASCs_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")

ASCs_D8v2<-FindMarkers(ASCs, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
ASCs_D8v2$"pct.1-pct.2" <- ASCs_D8v2$pct.1 - ASCs_D8v2$pct.2
ASCs_D8v2$"gene" <- rownames(ASCs_D8v2)
ASCs_D8v2$Celltype <- "ASCs"
ASCs_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(ASCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers.csv")
#save as txt file
write.table(ASCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(ASCs_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(ASCs_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
ASCs_D8v0<-FindMarkers(ASCs, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
ASCs_D8v0$"pct.1-pct.2" <- ASCs_D8v0$pct.1 - ASCs_D8v0$pct.2
ASCs_D8v0$"gene" <- rownames(ASCs_D8v0)
ASCs_D8v0$Celltype <- "ASCs"
ASCs_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(ASCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers.csv")
#save as txt file
write.table(ASCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
ASCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(ASCs_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(ASCs_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)



## pDCs
print("running analysis on pDCs")
pDCs<-subset(sal, idents = c("pDCs"))
#Get names of all genes in pDCs 
pDCs_bkg<-pDCs[["RNA"]]$counts
#remove non expressed genes
pDCs_bkg<-pDCs_bkg[rowSums(pDCs_bkg) > 0, ]
pDCs_bkg<-as.data.frame(pDCs_bkg)
pDCs_bkg<-rownames(pDCs)
pDCs_bkg<-as.data.frame(pDCs_bkg)
colnames(pDCs_bkg)<-"Gene"
#Add human geneID to the pDCs_bkg
pDCs_bkg<-left_join(pDCs_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
pDCs_bkg<-pDCs_bkg[!duplicated(pDCs_bkg$Gene), ]
write.table(pDCs_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(pDCs)<-"time_point"
print("table(Idents(pDCs))")
table(Idents(pDCs))
#Run FindMarkers on pDCs
pDCs_D2v0<-FindMarkers(object = pDCs, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
pDCs_D2v0$"pct.1-pct.2" <- pDCs_D2v0$pct.1 - pDCs_D2v0$pct.2
pDCs_D2v0$"gene" <- rownames(pDCs_D2v0)
pDCs_D2v0$Celltype <- "pDCs"
pDCs_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(pDCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers.csv")
#save as txt file
write.table(pDCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
pDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(pDCs_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(pDCs_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")

pDCs_D8v2<-FindMarkers(pDCs, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
pDCs_D8v2$"pct.1-pct.2" <- pDCs_D8v2$pct.1 - pDCs_D8v2$pct.2
pDCs_D8v2$"gene" <- rownames(pDCs_D8v2)
pDCs_D8v2$Celltype <- "pDCs"
pDCs_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(pDCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers.csv")
#save as txt file
write.table(pDCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(pDCs_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(pDCs_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
pDCs_D8v0<-FindMarkers(pDCs, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
pDCs_D8v0$"pct.1-pct.2" <- pDCs_D8v0$pct.1 - pDCs_D8v0$pct.2
pDCs_D8v0$"gene" <- rownames(pDCs_D8v0)
pDCs_D8v0$Celltype <- "pDCs"
pDCs_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(pDCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers.csv")
#save as txt file
write.table(pDCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
pDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(pDCs_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(pDCs_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)

## cDCs
print("running analysis on cDCs")
cDCs<-subset(sal, idents = c("cDCs"))
#Get names of all genes in cDCs 
cDCs_bkg<-cDCs[["RNA"]]$counts
#remove non expressed genes
cDCs_bkg<-cDCs_bkg[rowSums(cDCs_bkg) > 0, ]
cDCs_bkg<-as.data.frame(cDCs_bkg)
cDCs_bkg<-rownames(cDCs)
cDCs_bkg<-as.data.frame(cDCs_bkg)
colnames(cDCs_bkg)<-"Gene"
#Add human geneID to the cDCs_bkg
cDCs_bkg<-left_join(cDCs_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
cDCs_bkg<-cDCs_bkg[!duplicated(cDCs_bkg$Gene), ]
write.table(cDCs_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(cDCs)<-"time_point"
print("table(Idents(cDCs))")
table(Idents(cDCs))
#Run FindMarkers on cDCs
cDCs_D2v0<-FindMarkers(object = cDCs, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
cDCs_D2v0$"pct.1-pct.2" <- cDCs_D2v0$pct.1 - cDCs_D2v0$pct.2
cDCs_D2v0$"gene" <- rownames(cDCs_D2v0)
cDCs_D2v0$Celltype <- "cDCs"
cDCs_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(cDCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers.csv")
#save as txt file
write.table(cDCs_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
cDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(cDCs_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(cDCs_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
print("top50_genes_df2")
top50$"gene"<-as.character(top50$"gene")
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
print("sig2<-left_join(sig,ORG2.subset,by=c(gene=Gene))")
sig$"gene"<-as.character(sig$"gene")
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
print("pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c(gene=Gene))")
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
print("sig3<-left_join(sig2,pigrnaatlas2,by=c(gene=Gene))")
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")

cDCs_D8v2<-FindMarkers(cDCs, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
cDCs_D8v2$"pct.1-pct.2" <- cDCs_D8v2$pct.1 - cDCs_D8v2$pct.2
cDCs_D8v2$"gene" <- rownames(cDCs_D8v2)
cDCs_D8v2$Celltype <- "cDCs"
cDCs_D8v2$Comparison <- "8DPI vs 2DPI"
write.csv(cDCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers.csv")
#save as txt file
write.table(cDCs_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(cDCs_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(cDCs_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
cDCs_D8v0<-FindMarkers(cDCs, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
cDCs_D8v0$"pct.1-pct.2" <- cDCs_D8v0$pct.1 - cDCs_D8v0$pct.2
cDCs_D8v0$"gene" <- rownames(cDCs_D8v0)
cDCs_D8v0$Celltype <- "cDCs"
cDCs_D8v0$Comparison <- "8DPI vs 0DPI"
write.csv(cDCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers.csv")
#save as txt file
write.table(cDCs_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
cDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(cDCs_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(cDCs_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
if (nrow(sig) > 0) {
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))}
if (exists("sig2")) {
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) }

rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file (D2v0 no sig genes, so only D8v0 and D8v2)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D2v0_sig3,D8v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)



## CD4+ ab T cells
print("running analysis on CD4+ ab T cells")
Idents(sal)<-"celltypes_subclustres02"
CD4Tcells<-subset(sal, idents = c("CD4+ ab T cells"))
#Get names of all genes in CD4Tcells 
CD4Tcells_bkg<-CD4Tcells[["RNA"]]$counts
#remove non expressed genes
CD4Tcells_bkg<-CD4Tcells_bkg[rowSums(CD4Tcells_bkg) > 0, ]
CD4Tcells_bkg<-as.data.frame(CD4Tcells_bkg)
CD4Tcells_bkg<-rownames(CD4Tcells)
CD4Tcells_bkg<-as.data.frame(CD4Tcells_bkg)
colnames(CD4Tcells_bkg)<-"Gene"
#Add human geneID to the CD4Tcells_bkg
CD4Tcells_bkg<-left_join(CD4Tcells_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
CD4Tcells_bkg<-CD4Tcells_bkg[!duplicated(CD4Tcells_bkg$Gene), ]
write.table(CD4Tcells_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
Idents(CD4Tcells)<-"time_point"
print("table(Idents(CD4Tcells))")
table(Idents(CD4Tcells))
#Run FindMarkers on CD4Tcells
CD4Tcells_D2v0<-FindMarkers(object = CD4Tcells, slot="data", ident.1 ="D2", ident.2 ="D0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD4Tcells_D2v0$"pct.1-pct.2" <- CD4Tcells_D2v0$pct.1 - CD4Tcells_D2v0$pct.2
CD4Tcells_D2v0$"gene" <- rownames(CD4Tcells_D2v0)
CD4Tcells_D2v0$Celltype <- "CD4+ ab T cells"
CD4Tcells_D2v0$Comparison <- "2DPI vs 0DPI"
write.csv(CD4Tcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers.csv")
#save as txt file
write.table(CD4Tcells_D2v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers.txt", sep="\t",row.names=FALSE)
CD4Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(CD4Tcells_D2v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD4Tcells_D2v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 2####
print("Day 8 vs Day 2")

CD4Tcells_D8v2<-FindMarkers(CD4Tcells, slot="data",ident.1 ="D8", ident.2 ="D2" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD4Tcells_D8v2$"pct.1-pct.2" <- CD4Tcells_D8v2$pct.1 - CD4Tcells_D8v2$pct.2
CD4Tcells_D8v2$"gene" <- rownames(CD4Tcells_D8v2)
CD4Tcells_D8v2$"Celltype" <- "CD4+ ab T cells"
CD4Tcells_D8v2$"Comparison" <- "8DPI vs 2DPI"
write.csv(CD4Tcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers.csv")
#save as txt file
write.table(CD4Tcells_D8v2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers.txt", sep="\t",row.names=FALSE)
sig_pos <- subset(CD4Tcells_D8v2, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD4Tcells_D8v2, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)




###Day 8 vs Day 0####
print("Day 8 vs Day 0")
CD4Tcells_D8v0<-FindMarkers(CD4Tcells, slot="data",ident.1 ="D8",ident.2 ="D0",test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
CD4Tcells_D8v0$"pct.1-pct.2" <- CD4Tcells_D8v0$pct.1 - CD4Tcells_D8v0$pct.2
CD4Tcells_D8v0$"gene" <- rownames(CD4Tcells_D8v0)
CD4Tcells_D8v0$"Celltype" <- "CD4+ ab T cells"
CD4Tcells_D8v0$"Comparison" <- "8DPI vs 0DPI"
write.csv(CD4Tcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers.csv")
#save as txt file
write.table(CD4Tcells_D8v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers.txt", sep="\t",row.names=FALSE)
CD4Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(CD4Tcells_D8v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(CD4Tcells_D8v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
#Combine the positive and negative sig genes
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)


sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50$"gene"<-as.character(top50$"gene")
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig$"gene"<-as.character(sig$"gene")
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)

##merge all sig genes into one file
D8v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D2v0_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
D8v2_sig3<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
all_sig3<-rbind(D8v0_sig3,D2v0_sig3,D8v2_sig3)
write.table(all_sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_allComparisons_sig_genes.txt", sep="\t",row.names=FALSE)



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

cDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
if (exists("cDCs_D2v0") && nrow(cDCs_D2v0) > 0) {
cDCs_D2v0_2<-left_join(cDCs_D2v0,ORG2.subset,by=c("gene"="Gene"))
}
if (exists("cDCs_D2v0_2") && nrow(cDCs_D2v0_2) > 0) {
write.table(cDCs_D2v0_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE) 
}

cDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG.txt", sep = "\t", header = TRUE)
if (exists("cDCs_D8v2")) {
cDCs_D8v2_2<-left_join(cDCs_D8v2,ORG2.subset,by=c("gene"="Gene"))
}
if (exists("cDCs_D8v2_2 ") && nrow(cDCs_D8v2_2) > 0) {
write.table(cDCs_D8v2_2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)
}

#combine
if (exists("cDCs_D8v2_2") && exists("cDCs_D2v0_2")) {
cDCs_day <- rbind(cDCs_D8v2_2, cDCs_D2v0_2)
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

######################### Merge all DE gene files for paper #########################
Monocytes_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

Monocytes_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

Monocytes_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Monocytes_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

NKcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

NKcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

NKcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/NKcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

Bcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

Bcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

Bcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Bcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D2v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D8v2 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

CD2neg_GD_Tcells_D8v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD2neg_GD_Tcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

ASCs_D2v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

ASCs_D8v2 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

ASCs_D8v0 <-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/ASCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

pDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

pDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

pDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/pDCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

cDCs_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)


cDCs_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

cDCs_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/cDCs_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

CD4Tcells_D2v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D2v0_FindMarkers.txt", sep = "\t", header = TRUE)

CD4Tcells_D8v2<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v2_FindMarkers.txt", sep = "\t", header = TRUE)

CD4Tcells_D8v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/CD4Tcells_D8v0_FindMarkers.txt", sep = "\t", header = TRUE)

# Combine all data frames into one
all <- rbind(ASCs_D2v0, ASCs_D8v0, ASCs_D8v2,Bcells_D2v0, Bcells_D8v0, Bcells_D8v2,CD2neg_GD_Tcells_D2v0, CD2neg_GD_Tcells_D8v0, CD2neg_GD_Tcells_D8v2,CD4Tcells_D2v0, CD4Tcells_D8v0, CD4Tcells_D8v2,cDCs_D2v0, cDCs_D8v0, cDCs_D8v2,Monocytes_D2v0, Monocytes_D8v0, Monocytes_D8v2,NKcells_D2v0, NKcells_D8v0, NKcells_D8v2,pDCs_D2v0, pDCs_D8v0, pDCs_D8v2)
# Write the summary to a file
write.table(all, "/scRNAseq/Sal_5pigs_2026/DE/between_timepoints/Sal_FindMarkers_pairwise_DEGs_allCelltypes.txt", sep = "\t", quote = FALSE, row.names = TRUE)


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
#   [1] RcppAnnoy_0.0.22            splines_4.3.3              
#   [3] later_1.4.2                 ggplotify_0.1.2            
#   [5] bitops_1.0-9                polyclip_1.10-7            
#   [7] janitor_2.2.1               fastDummies_1.7.5          
#   [9] lifecycle_1.0.4             globals_0.16.3             
#  [11] lattice_0.22-6              MASS_7.3-60                
#  [13] MAST_1.28.0                 magrittr_2.0.3             
#  [15] plotly_4.10.4               httpuv_1.6.15              
#  [17] sctransform_0.4.1           spam_2.11-1                
#  [19] spatstat.sparse_3.1-0       reticulate_1.41.0.1        
#  [21] cowplot_1.1.3               pbapply_1.7-2              
#  [23] DBI_1.2.3                   abind_1.4-8                
#  [25] zlibbioc_1.48.2             GenomicRanges_1.54.1       
#  [27] Rtsne_0.17                  ggraph_2.2.1               
#  [29] RCurl_1.98-1.17             yulab.utils_0.2.0          
#  [31] tweenr_2.0.3                circlize_0.4.16            
#  [33] GenomeInfoDbData_1.2.11     ggrepel_0.9.6              
#  [35] irlba_2.3.5.1               listenv_0.9.1              
#  [37] spatstat.utils_3.1-2        tidytree_0.4.6             
#  [39] goftest_1.2-3               RSpectra_0.16-2            
#  [41] spatstat.random_3.3-2       fitdistrplus_1.2-2         
#  [43] parallelly_1.42.0           DelayedArray_0.28.0        
#  [45] codetools_0.2-20            ggforce_0.4.2              
#  [47] DOSE_3.28.2                 tidyselect_1.2.1           
#  [49] shape_1.4.6.1               aplot_0.2.5                
#  [51] farver_2.1.2                viridis_0.6.5              
#  [53] matrixStats_1.5.0           spatstat.explore_3.3-4     
#  [55] jsonlite_2.0.0              tidygraph_1.3.1            
#  [57] progressr_0.15.1            ggridges_0.5.6             
#  [59] survival_3.8-3              progress_1.2.3             
#  [61] tools_4.3.3                 treeio_1.26.0              
#  [63] ica_1.0-3                   Rcpp_1.0.14                
#  [65] glue_1.8.0                  SparseArray_1.2.4          
#  [67] gridExtra_2.3               MatrixGenerics_1.14.0      
#  [69] qvalue_2.34.0               GenomeInfoDb_1.38.8        
#  [71] withr_3.0.2                 fastmap_1.2.0              
#  [73] digest_0.6.37               gridGraphics_0.5-1         
#  [75] timechange_0.3.0            R6_2.6.1                   
#  [77] mime_0.13                   ggprism_1.0.5              
#  [79] colorspace_2.1-1            scattermore_1.2            
#  [81] GO.db_3.18.0                tensor_1.5                 
#  [83] spatstat.data_3.1-4         RSQLite_2.3.9              
#  [85] generics_0.1.3              data.table_1.17.0          
#  [87] prettyunits_1.2.0           S4Arrays_1.2.1             
#  [89] graphlayouts_1.2.2          httr_1.4.7                 
#  [91] htmlwidgets_1.6.4           scatterpie_0.2.4           
#  [93] uwot_0.2.3                  pkgconfig_2.0.3            
#  [95] gtable_0.3.6                blob_1.2.4                 
#  [97] lmtest_0.9-40               SingleCellExperiment_1.24.0
#  [99] XVector_0.42.0              shadowtext_0.1.4           
# [101] htmltools_0.5.8.1           dotCall64_1.2              
# [103] fgsea_1.28.0                scales_1.4.0               
# [105] png_0.1-8                   spatstat.univar_3.1-2      
# [107] snakecase_0.11.1            ggfun_0.1.8                
# [109] tzdb_0.4.0                  reshape2_1.4.4             
# [111] nlme_3.1-167                zoo_1.8-13                 
# [113] cachem_1.1.0                GlobalOptions_0.1.2        
# [115] KernSmooth_2.23-26          parallel_4.3.3             
# [117] miniUI_0.1.1.1              vipor_0.4.7                
# [119] HDO.db_0.99.1               ggrastr_1.0.2              
# [121] pillar_1.10.2               grid_4.3.3                 
# [123] vctrs_0.6.5                 RANN_2.6.2                 
# [125] promises_1.3.2              xtable_1.8-4               
# [127] cluster_2.1.8.1             beeswarm_0.4.0             
# [129] paletteer_1.6.0             cli_3.6.5                  
# [131] compiler_4.3.3              rlang_1.1.6                
# [133] crayon_1.5.3                future.apply_1.11.3        
# [135] rematch2_2.1.2              plyr_1.8.9                 
# [137] fs_1.6.6                    ggbeeswarm_0.7.2           
# [139] stringi_1.8.7               viridisLite_0.4.2          
# [141] deldir_2.0-4                BiocParallel_1.36.0        
# [143] babelgene_22.9              Biostrings_2.70.3          
# [145] lazyeval_0.2.2              spatstat.geom_3.3-5        
# [147] GOSemSim_2.28.1             RcppHNSW_0.6.0             
# [149] hms_1.1.3                   patchwork_1.3.0            
# [151] bit64_4.6.0-1               future_1.34.0              
# [153] KEGGREST_1.42.0             shiny_1.10.0               
# [155] SummarizedExperiment_1.32.0 ROCR_1.0-11                
# [157] igraph_2.1.4                memoise_2.0.1              
# [159] ggtree_3.10.1               fastmatch_1.1-6            
# [161] bit_4.6.0                   gson_0.1.0                 
# [163] ape_5.8-1                  