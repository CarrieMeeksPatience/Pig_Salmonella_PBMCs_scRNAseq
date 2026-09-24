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
Monocytes_bkg<-rownames(Monocytes_bkg)
Monocytes_bkg<-as.data.frame(Monocytes_bkg)
colnames(Monocytes_bkg)<-"Gene"
#Add human geneID to the Monocytes_bkg
Monocytes_bkg<-left_join(Monocytes_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
Monocytes_bkg<-Monocytes_bkg[!duplicated(Monocytes_bkg$Gene), ]
write.table(Monocytes_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

Idents(Monocytes)<-"sub.cluster.clustAll_res.0.09"
print("table(Idents(Monocytes))")
table(Idents(Monocytes))
#Run FindMarkers on Monocytes
Monocytes_C1v0<-FindMarkers(object = Monocytes, slot="data", ident.1 ="0_1", ident.2 ="0_0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Monocytes_C1v0$"pct.1-pct.2" <- Monocytes_C1v0$pct.1 - Monocytes_C1v0$pct.2
Monocytes_C1v0$"gene" <- rownames(Monocytes_C1v0)
Monocytes_C1v0$Celltype <- "Monocytes"
Monocytes_C1v0$Comparison <- "0_1 vs 0_0"

write.csv(Monocytes_C1v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers.csv")
#save as txt file
write.table(Monocytes_C1v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers.txt", sep="\t",row.names=FALSE)
Monocytes_C1v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Monocytes_C1v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Monocytes_C1v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_1v0_0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)


Monocytes_C0v1<-FindMarkers(object = Monocytes, slot="data", ident.1 ="0_0", ident.2 ="0_1" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
Monocytes_C0v1$"pct.1-pct.2" <- Monocytes_C0v1$pct.1 - Monocytes_C0v1$pct.2
Monocytes_C0v1$"gene" <- rownames(Monocytes_C0v1)
Monocytes_C0v1$Celltype <- "Monocytes"
Monocytes_C0v1$Comparison <- "0_0 vs 0_1"

write.csv(Monocytes_C0v1, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers.csv")
#save as txt file
write.table(Monocytes_C0v1, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers.txt", sep="\t",row.names=FALSE)
Monocytes_C0v1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(Monocytes_C0v1, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(Monocytes_C0v1, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_0_0v0_1_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)
Monocytes_all <- rbind(Monocytes_C1v0,Monocytes_C0v1)
table(Monocytes_all$Comparison)
#replace "0_1 vs 0_0" with "Monocyte 1" and "0_0 vs 0_1" with "Monocyte 0"
Monocytes_all$Comparison <- gsub("0_1 vs 0_0", "Monocyte 1", Monocytes_all$Comparison)
Monocytes_all$Comparison <- gsub("0_0 vs 0_1", "Monocyte 0", Monocytes_all$Comparison)
write.table(Monocytes_all, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/Monocytes_subcluster_pairwise_FindMarkers_all.txt", sep="\t",row.names=FALSE)

Idents(Monocytes)<-"sub.cluster.clustAll_res.0.1"
print("table(Idents(Monocytes))")
table(Idents(Monocytes))
base_dir <- "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/"
file_prefix <- "Monocytes_subRes01_subcluster_"

idents <- levels(Idents(Monocytes))
pairs <- combn(idents, 2, simplify = FALSE)

for(p in pairs) {
  a <- p[1]; b <- p[2]
  comp_tag <- paste0("Clusters", a, "v", b)            # e.g. Clusters0_0v0_1
  fname_base <- paste0(file_prefix, comp_tag)
  message("Running FindMarkers for ", a, " vs ", b, " -> ", comp_tag)
  mm <- tryCatch({
    FindMarkers(object = Monocytes, slot = "data", ident.1 = a, ident.2 = b, test.use = "MAST")
  }, error = function(e) {
    message("FindMarkers failed for ", comp_tag, ": ", e$message)
    return(NULL)
  })
  if(is.null(mm) || nrow(mm) == 0) next

  mm$`pct.1-pct.2` <- mm$pct.1 - mm$pct.2
  mm$gene <- rownames(mm)
  mm$ClusterComparison <- comp_tag

  out_all <- file.path(base_dir, paste0(fname_base, "_FindMarkers.txt"))
  write.table(mm, file = out_all, sep = "\t", row.names = FALSE)

  
  # cleanup iteration vars (optional)
  rm(mm)
}


idents <- levels(Idents(Monocytes))
pairs <- combn(idents, 2, simplify = FALSE)
print("Running FindMarkers for all cluster pairs in reverse order using for loop")
for(p in pairs) {
  a <- p[1]; b <- p[2]
  comp_tag <- paste0("Clusters", b, "v", a)            # e.g. Clusters0_0v0_1
  fname_base <- paste0(file_prefix, comp_tag)
  message("Running FindMarkers for ", b, " vs ", a, " -> ", comp_tag)
  mm <- tryCatch({
    FindMarkers(object = Monocytes, slot = "data", ident.1 = b, ident.2 = a, test.use = "MAST")
  }, error = function(e) {
    message("FindMarkers failed for ", comp_tag, ": ", e$message)
    return(NULL)
  })
  if(is.null(mm) || nrow(mm) == 0) next

  mm$`pct.1-pct.2` <- mm$pct.1 - mm$pct.2
  mm$gene <- rownames(mm)
  mm$ClusterComparison <- comp_tag

  out_all <- file.path(base_dir, paste0(fname_base, "_FindMarkers.txt"))
  write.table(mm, file = out_all, sep = "\t", row.names = FALSE)

  
  # cleanup iteration vars (optional)
  rm(mm)
}
# now combine the significant genes from all comparisons into one file for downstream analysis
files <- list.files(base_dir, pattern = "^Monocytes_subRes01_subcluster_.*_FindMarkers\\.txt$", full.names = TRUE)
# Combine all files into one data frame
all_sig <- do.call(rbind, lapply(files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t")
  return(df)
}))
write.table(all_sig, file = file.path(base_dir, "Monocytes_subRes01_subcluster_ALL_FindMarkers.txt"), sep = "\t", row.names = FALSE)


NKcells<-subset(sal, idents = c("NK cells"))
#Get names of all genes in NKcells 
NKcells_bkg<-NKcells[["RNA"]]$counts
#remove non expressed genes
NKcells_bkg<-NKcells_bkg[rowSums(NKcells_bkg) > 0, ]
NKcells_bkg<-as.data.frame(NKcells_bkg)
NKcells_bkg<-rownames(NKcells_bkg)
NKcells_bkg<-as.data.frame(NKcells_bkg)
colnames(NKcells_bkg)<-"Gene"
#Add human geneID to the NKcells_bkg
NKcells_bkg<-left_join(NKcells_bkg,ORG2.subset,by=c("Gene"="Gene"))
#remove duplicates
NKcells_bkg<-NKcells_bkg[!duplicated(NKcells_bkg$Gene), ]
write.table(NKcells_bkg,file="/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_bkg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

Idents(NKcells)<-"sub.cluster.clustAll_res.0.09"
print("table(Idents(NKcells))")
table(Idents(NKcells))
#Run FindMarkers on NKcells
NKcells_C1v0<-FindMarkers(object = NKcells, slot="data", ident.1 ="4_1", ident.2 ="4_0" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
NKcells_C1v0$"pct.1-pct.2" <- NKcells_C1v0$pct.1 - NKcells_C1v0$pct.2
NKcells_C1v0$"gene" <- rownames(NKcells_C1v0)
NKcells_C1v0$Celltype <- "NKcells"
NKcells_C1v0$Comparison <- "4_1 vs 4_0"

write.csv(NKcells_C1v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers.csv")
#save as txt file
write.table(NKcells_C1v0, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers.txt", sep="\t",row.names=FALSE)
NKcells_C1v0<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(NKcells_C1v0, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(NKcells_C1v0, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_1v4_0_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)


NKcells_C0v1<-FindMarkers(object = NKcells, slot="data", ident.1 ="4_0", ident.2 ="4_1" ,test.use="MAST")
#add column "pct.1-pct.2" to the data frame, if the value is positive, the gene is more highly expressed in time_point 1, if the value is negative, the gene is more highly expressed in time_point 2
NKcells_C0v1$"pct.1-pct.2" <- NKcells_C0v1$pct.1 - NKcells_C0v1$pct.2
NKcells_C0v1$"gene" <- rownames(NKcells_C0v1)
NKcells_C0v1$Celltype <- "NKcells"
NKcells_C0v1$Comparison <- "4_0 vs 4_1"

write.csv(NKcells_C0v1, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers.csv")
#save as txt file
write.table(NKcells_C0v1, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers.txt", sep="\t",row.names=FALSE)
NKcells_C0v1<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers.txt", sep = "\t", header = TRUE)
sig_pos <- subset(NKcells_C0v1, (avg_log2FC >= 0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig_neg <- subset(NKcells_C0v1, (avg_log2FC <= -0.25 & p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
sig <- rbind(sig_pos, sig_neg)
write.table(sig, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG.txt", sep="\t",row.names=FALSE)

sig<-read.table("/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG.txt", sep = "\t", header = TRUE)

#Select the top 50 genes(sig p.adjust & highest Log2FC)
top50<-sig %>%top_n(50,avg_log2FC)
write.table(top50, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG_top50_log2FC_genes_2026_04_16.txt", sep="\t",row.names=FALSE)

#Add pig & human esembl ids to the top50 genes
top50_genes_df2<-left_join(top50,ORG2.subset,by=c("gene"="Gene"))
sig2<-left_join(sig,ORG2.subset,by=c("gene"="Gene"))
write.table(sig2, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG_duplicateGenes_2026_04_16.txt", sep="\t",row.names=FALSE)

#remove rows with duplicated genes
top50_genes_df2<-top50_genes_df2[!duplicated(top50_genes_df2$Gene.stable.ID),]
sig2<-sig2[!duplicated(sig2$Gene.stable.ID),]
#Merge top50_genes_df with pigrnaatlas2
pigrnaatlas3<-left_join(top50_genes_df2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(pigrnaatlas3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG_top50genes_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
sig3<-left_join(sig2,pigrnaatlas2,by=c("gene"="Gene"))
write.table(sig3, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_4_0v4_1_FindMarkers_SIG_pigrnaatlas_2026_04_16.txt", sep="\t",row.names=FALSE)
rm(sig,sig2,sig3,top50,top50_genes_df2,pigrnaatlas3)
NKcells_all <- rbind(NKcells_C1v0,NKcells_C0v1)
table(NKcells_all$Comparison)
#replace "4_1 vs 4_0" with "NKcell 1" and "4_0 vs 4_1" with "NKcell 0"
NKcells_all$Comparison <- gsub("4_1 vs 4_0", "NKcell 1", NKcells_all$Comparison)
NKcells_all$Comparison <- gsub("4_0 vs 4_1", "NKcell 0", NKcells_all$Comparison)
write.table(NKcells_all, file = "/scRNAseq/Sal_5pigs_2026/DE/between_subcluster/NKcells_subcluster_pairwise_FindMarkers_all.txt", sep="\t",row.names=FALSE)

