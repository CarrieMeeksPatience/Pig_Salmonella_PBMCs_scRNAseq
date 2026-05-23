.libPaths()
#[1] "/micromamba/envs/seurat+milo/lib/R/library"

{set.seed(123)
library(Seurat)
library(ggplot2)
library(dplyr)
library(STACAS)
library(scIntegrationMetrics)
library(tidyr)
library(scCustomize)
library(clustree)
library(RColorBrewer)
library(PCAtools)
library(findPC)}

# Read in seurat object
All<- readRDS("/scRNAseq/Sal_5pigs_2026/Normalization/5pigs_salmonella_Filtered_ScranNorm_Regressed_PCA_2026_04_15.rds")
## 1. Retrieve full list of signatures for human
#hs.sign <- GetSignature(SignatuR$Hs)
```
#carrie's version 
##annotation<-/PIPseq_PIPseeker/ALL/gtf/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.csv
##gene_block_CSM<-filter(annotation,gene_biotype ==c("pseudogene") )
##gene_block_CSM<-gene_block_CSM$gene_name
##s.genes <- cc.genes$s.genes
##g2m.genes <- cc.genes$g2m.genes

## 2. Then we tell STACAS to exclude specific genes from HVG used in integration, using the parameter `genesBlockList`
#my.genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),GetSignature(SignatuR$Hs$Compartments))

#Carrie's version
##gene_block_list_CSM<-list("pseudogene" = gene_block_CSM,"s.genes"= s.genes, "g2m.genes"=g2m.genes )
##save<-dput( gene_block_list_CSM, file="/scRNAseq/All_Flu_902/Integration_test/gene_block_list_CSM_2024_06_17.R")
gene_block_list_CSM<-dget("/scRNAseq/All_Flu_902/Integration_test/gene_block_list_CSM_2024_06_17.R")
## Try integration via Animal PC1:11 ##
# subset whole dataset 

All1<-All
ndim=14
nfeatures<-VariableFeatures(All)
All1 <- RunUMAP(All1,dims=1:ndim)
colorblind_palette <- ColorBlind_Pal()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d14_1280HVF.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All1, group.by = "SampleID") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d14_1280HVF")
A2<-DimPlot_scCustom(All1, group.by = "Animal") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d14_1280HVF")
A1|A2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d14_1280HVF_timepoint.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All1, group.by = "Animal", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d14_1280HVF")
A2<-DimPlot_scCustom(All1, group.by = "time_point", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d14_1280HVF")
A1|A2
dev.off()

## 1. Find integration anchors between datasets/batches
obj.list1 <- SplitObject(All1, split.by = "Animal")
stacas_anchors1 <- FindAnchors.STACAS(obj.list1, genesBlockList=gene_block_list_CSM,anchor.features = nfeatures,dims = 1:ndim)

## 2. Guide tree for integration order
st1 <- SampleTree.STACAS(anchorset = stacas_anchors1,obj.names = names(obj.list1))

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/sample_tree_Animal_ScranNorm_Regressed_PCA_d14_1280HVF.pdf")
SampleTree.STACAS(anchorset = stacas_anchors1,obj.names = names(obj.list1))
dev.off()
# Are samples from the same sequencing technology or from similar tissues clustering together? Different hierarchical clustering methods are available as `hclust.methods` parameter.

## 3. Dataset integration
All1_integrated <- IntegrateData.STACAS(stacas_anchors1,sample.tree = st1,dims=1:ndim)
  #Warning: Layer counts isn't present in the assay object; returning NULL
  #Warning: Layer counts isn't present in the assay object; returning NULL
All1_s<-saveRDS(All1_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF.rds")
#Calculate low-dimensional embeddings and visualize integration results in UMAP
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# get cell cycle scores
All1_integrated <- CellCycleScoring(All1_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regressed out cycling genes
All1_integrated  <- ScaleData(All1_integrated , vars.to.Regressed = c("S.Score", "G2M.Score"), features = rownames(All1_integrated))
All1_integrated<- RunPCA(All1_integrated, npcs = 100, verbose = TRUE)
#Elbow plot
pdf("/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF_elbow.pdf")
ElbowPlot(All1_integrated, ndims = 40) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d14_1280HVF")
dev.off()
#Run UMAP
All1_integrated <- RunUMAP(All1_integrated, dims = 1:14)
All1_s<-saveRDS(All1_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF.rds")


pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF.pdf", width=25, height=15)
P1<-DimPlot_scCustom(All1_integrated, group.by = "SampleID") +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d14_1280HVF")
P2<-DimPlot_scCustom(All1_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d14_1280HVF.")
P1|P2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF_v2.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All1_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d14_1280HVF")
P2<-DimPlot_scCustom(All1_integrated, group.by = "time_point", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d14_1280HVF")
P1|P2
dev.off()
rm(ndim)


#Try with PC1:9
All2<-All
ndim=9
nfeatures<-VariableFeatures(All)
All2 <- RunUMAP(All2,dims=1:ndim)
colorblind_palette <- ColorBlind_Pal()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d9_1280HVF.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All2, group.by = "SampleID") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d9_1280HVF")
A2<-DimPlot_scCustom(All2, group.by = "Animal") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d9_1280HVF")
A1|A2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d9_1280HVF_timepoint.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All2, group.by = "Animal", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d9_1280HVF")
A2<-DimPlot_scCustom(All2, group.by = "time_point", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d9_1280HVF")
A1|A2
dev.off()

## 1. Find integration anchors between datasets/batches
obj.list2 <- SplitObject(All2, split.by = "Animal")
stacas_anchors2 <- FindAnchors.STACAS(obj.list2, genesBlockList=gene_block_list_CSM,anchor.features = nfeatures,dims = 1:ndim)

## 2. Guide tree for integration order
st2 <- SampleTree.STACAS(anchorset = stacas_anchors2,obj.names = names(obj.list2))

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/sample_tree_Animal_ScranNorm_Regressed_PCA_d9_1280HVF.pdf")
SampleTree.STACAS(anchorset = stacas_anchors2,obj.names = names(obj.list2))
dev.off()
# Are samples from the same sequencing technology or from similar tissues clustering together? Different hierarchical clustering methods are available as `hclust.methods` parameter.

## 3. Dataset integration
All2_integrated <- IntegrateData.STACAS(stacas_anchors2,sample.tree = st2,dims=1:ndim)
  #Warning: Layer counts isn't present in the assay object; returning NULL
  #Warning: Layer counts isn't present in the assay object; returning NULL
All2_s<-saveRDS(All2_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d9_1280HVF.rds")
#Calculate low-dimensional embeddings and visualize integration results in UMAP
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# get cell cycle scores
All2_integrated <- CellCycleScoring(All2_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regressed out cycling genes
All2_integrated  <- ScaleData(All2_integrated , vars.to.Regressed = c("S.Score", "G2M.Score"), features = rownames(All2_integrated))
All2_integrated<- RunPCA(All2_integrated, npcs = 100, verbose = TRUE)
#Elbow plot
pdf("/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d9_1280HVF_elbow.pdf")
ElbowPlot(All2_integrated, ndims = 40) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d9_1280HVF")
dev.off()
#Run UMAP
All2_integrated <- RunUMAP(All2_integrated, dims = 1:9)
All2_s<-saveRDS(All2_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d9_1280HVF.rds")


pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d9_1280HVF.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All2_integrated, group.by = "SampleID") +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d9_1280HVF")
P2<-DimPlot_scCustom(All2_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d9_1280HVF")
P1|P2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d9_1280HVF_v2.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All2_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d9_1280HVF")
P2<-DimPlot_scCustom(All2_integrated, group.by = "time_point", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d9_1280HVF")
P1|P2
dev.off()
rm(ndim)

#Try with PC1:10
All3<-All
ndim=10
nfeatures<-VariableFeatures(All)
All3 <- RunUMAP(All3,dims=1:ndim)
colorblind_palette <- ColorBlind_Pal()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d10_1280HVF.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All3, group.by = "SampleID") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d10_1280HVF")
A2<-DimPlot_scCustom(All3, group.by = "Animal") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d10_1280HVF")
A1|A2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d10_1280HVF_timepoint.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All3, group.by = "Animal", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d10_1280HVF")
A2<-DimPlot_scCustom(All3, group.by = "time_point", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d10_1280HVF")
A1|A2
dev.off()

## 1. Find integration anchors between datasets/batches
obj.list3 <- SplitObject(All3, split.by = "Animal")
stacas_anchors3 <- FindAnchors.STACAS(obj.list3, genesBlockList=gene_block_list_CSM,anchor.features = nfeatures,dims = 1:ndim)

## 2. Guide tree for integration order
st3 <- SampleTree.STACAS(anchorset = stacas_anchors3,obj.names = names(obj.list3))

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/sample_tree_Animal_ScranNorm_Regressed_PCA_d10_1280HVF.pdf")
SampleTree.STACAS(anchorset = stacas_anchors3,obj.names = names(obj.list3))
dev.off()
# Are samples from the same sequencing technology or from similar tissues clustering together? Different hierarchical clustering methods are available as `hclust.methods` parameter.

## 3. Dataset integration
All3_integrated <- IntegrateData.STACAS(stacas_anchors3,sample.tree = st3,dims=1:ndim)
  #Warning: Layer counts isn't present in the assay object; returning NULL
  #Warning: Layer counts isn't present in the assay object; returning NULL
All3_s<-saveRDS(All3_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF.rds")
#Calculate low-dimensional embeddings and visualize integration results in UMAP
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# get cell cycle scores
All3_integrated <- CellCycleScoring(All3_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regressed out cycling genes
All3_integrated  <- ScaleData(All3_integrated , vars.to.Regressed = c("S.Score", "G2M.Score"), features = rownames(All3_integrated))
All3_integrated<- RunPCA(All3_integrated, npcs = 100, verbose = TRUE)
#Elbow plot
pdf("/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF_elbow.pdf")
ElbowPlot(All3_integrated, ndims = 40) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1280HVF")
dev.off()
#Run UMAP
All3_integrated <- RunUMAP(All3_integrated, dims = 1:10)
All3_s<-saveRDS(All3_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF.rds")


pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All3_integrated, group.by = "SampleID") +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1280HVF")
P2<-DimPlot_scCustom(All3_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1280HVF.")
P1|P2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF_v2.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All3_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1280HVF")
P2<-DimPlot_scCustom(All3_integrated, group.by = "time_point", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1280HVF")
P1|P2
dev.off()
rm(ndim)

#Try with PC1:11
All4<-All
ndim=11
nfeatures<-VariableFeatures(All)
All4 <- RunUMAP(All4,dims=1:ndim)
colorblind_palette <- ColorBlind_Pal()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d11_1280HVF.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All4, group.by = "SampleID") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d11_1280HVF")
A2<-DimPlot_scCustom(All4, group.by = "Animal") + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d11_1280HVF")
A1|A2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/nonIntegrated_dimPlot_ScranNorm_Regressed_PCA_d11_1280HVF_timepoint.pdf", width = 25 , height=15)
A1<-DimPlot_scCustom(All4, group.by = "Animal", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d11_1280HVF")
A2<-DimPlot_scCustom(All4, group.by = "time_point", colors_use = colorblind_palette) + theme(aspect.ratio = 1) + ggtitle("Non-integrated ScranNorm_Regressed_PCA_d11_1280HVF")
A1|A2
dev.off()

## 1. Find integration anchors between datasets/batches
obj.list4 <- SplitObject(All4, split.by = "Animal")
stacas_anchors4 <- FindAnchors.STACAS(obj.list4, genesBlockList=gene_block_list_CSM,anchor.features = nfeatures,dims = 1:ndim)

## 2. Guide tree for integration order
st4 <- SampleTree.STACAS(anchorset = stacas_anchors4,obj.names = names(obj.list4))

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/sample_tree_Animal_ScranNorm_Regressed_PCA_d11_1280HVF.pdf")
SampleTree.STACAS(anchorset = stacas_anchors4,obj.names = names(obj.list4))
dev.off()
# Are samples from the same sequencing technology or from similar tissues clustering together? Different hierarchical clustering methods are available as `hclust.methods` parameter.

## 3. Dataset integration
All4_integrated <- IntegrateData.STACAS(stacas_anchors4,sample.tree = st4,dims=1:ndim)
  #Warning: Layer counts isn't present in the assay object; returning NULL
  #Warning: Layer counts isn't present in the assay object; returning NULL
All4_s<-saveRDS(All4_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d11_1280HVF.rds")
#Calculate low-dimensional embeddings and visualize integration results in UMAP
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# get cell cycle scores
All4_integrated <- CellCycleScoring(All4_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regressed out cycling genes
All4_integrated  <- ScaleData(All4_integrated , vars.to.Regressed = c("S.Score", "G2M.Score"), features = rownames(All4_integrated))
All4_integrated<- RunPCA(All4_integrated, npcs = 100, verbose = TRUE)
#Elbow plot
pdf("/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d11_1280HVF_elbow.pdf")
ElbowPlot(All4_integrated, ndims = 40) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d11_1280HVF")
dev.off()
#Run UMAP
All4_integrated <- RunUMAP(All4_integrated, dims = 1:11)
All4_s<-saveRDS(All4_integrated, file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d11_1280HVF.rds")


pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d11_1280HVF.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All4_integrated, group.by = "SampleID") +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d11_1280HVF")
P2<-DimPlot_scCustom(All4_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d11_1280HVF.")
P1|P2
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Integration/STACAS_UMAP_via_Animal_ScranNorm_Regressed_PCA_d11_1280HVF_v2.pdf", width=25, height =15)
P1<-DimPlot_scCustom(All4_integrated, group.by = "Animal", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d11_1280HVF")
P2<-DimPlot_scCustom(All4_integrated, group.by = "time_point", colors_use = colorblind_palette) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d11_1280HVF")
P1|P2
dev.off()
rm(ndim)

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
#  [1] tidyr_1.3.1                 scIntegrationMetrics_1.1
#  [3] STACAS_2.2.2                patchwork_1.3.0
#  [5] pheatmap_1.0.12             igraph_2.1.4
#  [7] bluster_1.12.0              findPC_1.0
#  [9] PCAtools_2.14.0             ggrepel_0.9.6
# [11] RColorBrewer_1.1-3          clustree_0.5.1
# [13] ggraph_2.2.1                scCustomize_3.0.1
# [15] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [17] Biobase_2.62.0              GenomicRanges_1.54.1
# [19] GenomeInfoDb_1.38.8         IRanges_2.36.0
# [21] S4Vectors_0.40.2            BiocGenerics_0.48.1
# [23] MatrixGenerics_1.14.0       matrixStats_1.5.0
# [25] dplyr_1.1.4                 ggplot2_3.5.2
# [27] Matrix_1.6-5                SeuratWrappers_0.3.5
# [29] pbmcsca.SeuratData_3.0.0    pbmcref.SeuratData_1.0.0
# [31] pbmc3k.SeuratData_3.1.4     ifnb.SeuratData_3.1.0
# [33] cbmc.SeuratData_3.1.4       SeuratData_0.2.2.9001
# [35] Seurat_5.2.1                SeuratObject_5.0.2
# [37] sp_2.2-0

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.3.3
#   [3] later_1.4.2               prismatic_1.1.2
#   [5] bitops_1.0-9              tibble_3.2.1
#   [7] R.oo_1.27.0               polyclip_1.10-7
#   [9] janitor_2.2.1             fastDummies_1.7.5
#  [11] lifecycle_1.0.4           globals_0.16.3
#  [13] lattice_0.22-6            MASS_7.3-60
#  [15] magrittr_2.0.3            plotly_4.10.4
#  [17] remotes_2.5.0             httpuv_1.6.15
#  [19] sctransform_0.4.1         spam_2.11-1
#  [21] spatstat.sparse_3.1-0     reticulate_1.41.0.1
#  [23] cowplot_1.1.3             pbapply_1.7-2
#  [25] lubridate_1.9.4           abind_1.4-8
#  [27] zlibbioc_1.48.2           Rtsne_0.17
#  [29] purrr_1.0.4               R.utils_2.13.0
#  [31] RCurl_1.98-1.17           tweenr_2.0.3
#  [33] rappdirs_0.3.3            circlize_0.4.16
#  [35] GenomeInfoDbData_1.2.11   irlba_2.3.5.1
#  [37] listenv_0.9.1             spatstat.utils_3.1-2
#  [39] vegan_2.6-10              goftest_1.2-3
#  [41] RSpectra_0.16-2           dqrng_0.4.1
#  [43] spatstat.random_3.3-2     fitdistrplus_1.2-2
#  [45] parallelly_1.42.0         permute_0.9-7
#  [47] DelayedMatrixStats_1.24.0 codetools_0.2-20
#  [49] DelayedArray_0.28.0       ggforce_0.4.2
#  [51] tidyselect_1.2.1          shape_1.4.6.1
#  [53] farver_2.1.2              ScaledMatrix_1.10.0
#  [55] viridis_0.6.5             spatstat.explore_3.3-4
#  [57] jsonlite_2.0.0            BiocNeighbors_1.20.2
#  [59] tidygraph_1.3.1           progressr_0.15.1
#  [61] ggridges_0.5.6            survival_3.8-3
#  [63] tools_4.3.3               ica_1.0-3
#  [65] Rcpp_1.0.14               glue_1.8.0
#  [67] gridExtra_2.3             SparseArray_1.2.4
#  [69] mgcv_1.9-1                withr_3.0.2
#  [71] BiocManager_1.30.25       fastmap_1.2.0
#  [73] digest_0.6.37             rsvd_1.0.5
#  [75] timechange_0.3.0          R6_2.6.1
#  [77] mime_0.13                 ggprism_1.0.5
#  [79] colorspace_2.1-1          scattermore_1.2
#  [81] tensor_1.5                spatstat.data_3.1-4
#  [83] R.methodsS3_1.8.2         generics_0.1.3
#  [85] data.table_1.17.0         graphlayouts_1.2.2
#  [87] httr_1.4.7                htmlwidgets_1.6.4
#  [89] S4Arrays_1.2.1            uwot_0.2.3
#  [91] pkgconfig_2.0.3           gtable_0.3.6
#  [93] lmtest_0.9-40             XVector_0.42.0
#  [95] htmltools_0.5.8.1         dotCall64_1.2
#  [97] scales_1.4.0              png_0.1-8
#  [99] spatstat.univar_3.1-2     snakecase_0.11.1
# [101] reshape2_1.4.4            nlme_3.1-167
# [103] cachem_1.1.0              zoo_1.8-13
# [105] GlobalOptions_0.1.2       stringr_1.5.1
# [107] KernSmooth_2.23-26        parallel_4.3.3
# [109] miniUI_0.1.1.1            vipor_0.4.7
# [111] ggrastr_1.0.2             pillar_1.10.2
# [113] grid_4.3.3                vctrs_0.6.5
# [115] RANN_2.6.2                promises_1.3.2
# [117] BiocSingular_1.18.0       beachmat_2.18.1
# [119] xtable_1.8-4              cluster_2.1.8.1
# [121] beeswarm_0.4.0            paletteer_1.6.0
# [123] cli_3.6.5                 compiler_4.3.3
# [125] rlang_1.1.6               crayon_1.5.3
# [127] future.apply_1.11.3       labeling_0.4.3
# [129] rematch2_2.1.2            plyr_1.8.9
# [131] forcats_1.0.0             ggbeeswarm_0.7.2
# [133] stringi_1.8.7             BiocParallel_1.36.0
# [135] viridisLite_0.4.2         deldir_2.0-4
# [137] lazyeval_0.2.2            spatstat.geom_3.3-5
# [139] RcppHNSW_0.6.0            sparseMatrixStats_1.14.0
# [141] future_1.34.0             shiny_1.10.0
# [143] ROCR_1.0-11               memoise_2.0.1
