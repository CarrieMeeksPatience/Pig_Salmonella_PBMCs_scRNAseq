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

Jayne_NK_markers<- c("PAX5","CD19","CD79A","TCF4","SLA-DRB1","SLA-DRA","CR2","IRF8","CD86","TYROBP","HCST","KLRB1", "KLRK1","PRF1","CD6","CD5","CD8A","CD8B","CD4","TRDC","CD2","CD3E","CD3G","PRDM1","IRF4","CD93","CLEC12A","XBP1","FCER1A","FLT3","TLR4","NLRP3","CSF1R","CD163","CD14","SIRPA", "IL2RB" , "FCGR3A", "NCAM1", "CD69", "KLRD1", "COX6A2", "KIR2DL4", "NKG7", "NCR1", "ZMAT4")


All1_integrated<-readRDS(file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d14_1280HVF.rds")
#Now Cluster


#Try PCs1:10
All3_integrated<-readRDS(file ="/scRNAseq/Sal_5pigs_2026/Integration/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_1280HVF.rds")
All3_integrated <- FindNeighbors(All3_integrated,dims = 1:10)

{All3_integrated <- FindClusters(All3_integrated, resolution = 0.1)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.2)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.3)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.4)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.5)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.6)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.7)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.8)
All3_integrated <- FindClusters(All3_integrated, resolution = 0.9)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.0)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.1)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.2)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.3)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.4)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.5)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.6)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.7)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.8)
All3_integrated <- FindClusters(All3_integrated, resolution = 1.9)
All3_integrated <- FindClusters(All3_integrated, resolution = 2.0)
All3_integrated <- FindClusters(All3_integrated, resolution = 2.1)
All3_integrated <- FindClusters(All3_integrated, resolution = 2.2)}

All3_integrated <- RunUMAP(All3_integrated, dims = 1:10)

saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_15.rds")

pdf(file="/scRNAseq/Sal_5pigs_2026//clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clusterTree.pdf", width=20, height=20)
clustree(All3_integrated,  prefix = "integrated_snn_res.")
dev.off()

DefaultAssay(All3_integrated) <- "RNA"
Idents(All3_integrated) <- "integrated_snn_res.0.1"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_res01_dotplot.pdf", width=6, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
dev.off()
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_res01_dotplot_horizontal.pdf", width=12, height=4)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,x_lab_rotate=45,flip = TRUE,grid_color = "grey")
dev.off()
Idents(All3_integrated) <- "integrated_snn_res.0.2"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_res02_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
dev.off()
Idents(All3_integrated) <- "integrated_snn_res.0.3"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_res03_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
dev.off()
saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_15.rds")

All3_integrated$"integrated_snn_res.0.6"<-NULL
All3_integrated$"integrated_snn_res.0.7"<-NULL
All3_integrated$"integrated_snn_res.0.8"<-NULL
All3_integrated$"integrated_snn_res.0.9"<-NULL
All3_integrated$"integrated_snn_res.1"<-NULL
All3_integrated$"integrated_snn_res.1.1"<-NULL
All3_integrated$"integrated_snn_res.1.2"<-NULL
All3_integrated$"integrated_snn_res.1.3"<-NULL
All3_integrated$"integrated_snn_res.1.4"<-NULL
All3_integrated$"integrated_snn_res.1.5"<-NULL
All3_integrated$"integrated_snn_res.1.6"<-NULL
All3_integrated$"integrated_snn_res.1.7"<-NULL
All3_integrated$"integrated_snn_res.1.8"<-NULL
All3_integrated$"integrated_snn_res.1.9"<-NULL
All3_integrated$"integrated_snn_res.2"<-NULL
All3_integrated$"integrated_snn_res.2.1"<-NULL
All3_integrated$"integrated_snn_res.2.2"<-NULL

pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clusterTree_short.pdf")
clustree(All3_integrated,  prefix = "integrated_snn_res.")
dev.off()


# I like PC1:10 the best, so now I will subcluster the PC1:10 object
All3_integrated <- readRDS("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_15.rds")
DefaultAssay(All3_integrated) <- "integrated"
Idents(All3_integrated) <- "integrated_snn_res.0.1"
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 0,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust0_res.0.2",
  resolution = 0.2,
  algorithm = 1
)

All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 1,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust1_res.0.2",
  resolution = 0.2,
  algorithm = 1
)

All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 2,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust2_res.0.2",
  resolution = 0.2,
  algorithm = 1
)

All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 3,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust3_res.0.2",
  resolution = 0.2,
  algorithm = 1
)

All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 4,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust4_res.0.2",
  resolution = 0.2,
  algorithm = 1
)

All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.05",
  resolution = 0.05,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.075",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.08",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.09",
  resolution = 0.075,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.1",
  resolution = 0.1,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.125",
  resolution = 0.125,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.15",
  resolution = 0.15,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.1525",
  resolution = 0.1525,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.155",
  resolution = 0.155,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.16",
  resolution = 0.16,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.175",
  resolution = 0.175,
  algorithm = 1
)
All3_integrated<- FindSubCluster(
  All3_integrated,
  cluster = 5,
  graph.name= "integrated_snn",
  subcluster.name = "sub.cluster.clust5_res.0.2",
  resolution = 0.2,
  algorithm = 1
)
#now combine the subcluster resolutions into one meta.data column for each resolution

All3_integrated$sub.cluster.clustAll_res.0.05 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.05")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.05[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.075 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.075")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.075[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.08 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.08")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.08[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.09 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.09")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.09[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.1 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.1")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.1[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.125 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.125")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.125[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.15 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.15")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.15[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.1525 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.1525")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.1525[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.155 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.155")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.155[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.16 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.16")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.16[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.175 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.175")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.175[idx] <- vals[idx]
}

All3_integrated$sub.cluster.clustAll_res.0.2 <- NA_character_
for(i in 0:5) {
  col_name <- paste0("sub.cluster.clust", i, "_res.0.2")
  if(! col_name %in% colnames(All3_integrated@meta.data)) {
    message("Skipping missing meta column: ", col_name)
    next
  }
  vals <- as.character(All3_integrated@meta.data[[col_name]])
  vals_unique <- unique(vals)
  vals_keep <- vals_unique[!is.na(vals_unique) & grepl("_", vals_unique)]
  if(length(vals_keep) == 0) next
  idx <- which(vals %in% vals_keep)
  All3_integrated$sub.cluster.clustAll_res.0.2[idx] <- vals[idx]
}
saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")
saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")

pdf("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclustered_clustAll_clustree.pdf", width=15,height=9)
clustree(All3_integrated,  prefix = "sub.cluster.clustAll_res.")
dev.off()

DefaultAssay(All3_integrated) <- "RNA"
Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.1"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres01_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.15"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres015_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.05"
Monocytes <- subset(All3_integrated, idents=c("0_0"))
 sal5_mono_gene <- c("HDAC9","ERBIN","AKT3","ANKRD17","NEK7","BIRC3","NFKBIA","IKBKB","DOCK4","ITGA4","APP","ITGAL","MERTK","CD74","S100A10","FCGR3A","TNFRSF1B","RACK1","BST2","B2M","CXCL10","GBP2","TGFB1","GRN","MYD88","STAT1","GBP1","IRF8","LYZ","STAT3","LTF","SOD2","CD163","MAPK14","S100A8","CD14","S100A9","FOS")
Idents(Monocytes) <- "sub.cluster.clustAll_res.0.09"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres009_mono_dotplot.pdf", width=5, height=12)
Clustered_DotPlot(Monocytes, sal5_mono_gene,grid_color = "grey")
dev.off()
Idents(Monocytes) <- "sub.cluster.clustAll_res.0.1"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres01_mono_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(Monocytes, sal5_mono_gene,grid_color = "grey")
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.05"
NKcells <- subset(All3_integrated, idents=c("4_0"))
 sal5_NK_gene <- c("FCER1G","B2M","CD2","CD74","CORO1A","HCST","HSPA8","KLRK1","NKG7","RAC2","S100A13","TYROBP","ADGRG1","CARD11","CBLB","IKZF3","IL12RB2","ITPKB","RORA","RUNX1","TOX","ZEB2","PPP3CA","DOCK8","PTK2B","STAT3","SUN2","TGFBR2","WNK1","AKT3","PTPRC","ITGA4")
Idents(NKcells) <- "sub.cluster.clustAll_res.0.09"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres009_nk_dotplot.pdf", width=5, height=12)
Clustered_DotPlot(NKcells, sal5_NK_gene,grid_color = "grey")
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.05"
Bcells <- subset(All3_integrated, idents=c("1_0","1_1"))
Idents(Bcells) <- "sub.cluster.clustAll_res.0.2"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres02_bcell_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(Bcells, Jayne_NK_markers,grid_color = "grey")
dev.off()
Idents(Bcells) <- "sub.cluster.clustAll_res.0.16"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres016_bcell_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(Bcells, Jayne_NK_markers,grid_color = "grey")
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.2"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclusterres02_dotplot.pdf", width=10, height=12)
Clustered_DotPlot(All3_integrated, Jayne_NK_markers,grid_color = "grey")
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
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] findPC_1.0               PCAtools_2.14.0          ggrepel_0.9.6
#  [4] RColorBrewer_1.1-3       clustree_0.5.1           ggraph_2.2.1
#  [7] scCustomize_3.0.1        tidyr_1.3.1              scIntegrationMetrics_1.1
# [10] STACAS_2.2.2             dplyr_1.1.4              ggplot2_3.5.2
# [13] Seurat_5.2.1             SeuratObject_5.0.2       sp_2.2-0

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.3.3
#   [3] later_1.4.2               prismatic_1.1.2
#   [5] tibble_3.2.1              R.oo_1.27.0
#   [7] polyclip_1.10-7           janitor_2.2.1
#   [9] fastDummies_1.7.5         lifecycle_1.0.4
#  [11] doParallel_1.0.17         globals_0.16.3
#  [13] lattice_0.22-6            MASS_7.3-60
#  [15] backports_1.5.0           magrittr_2.0.3
#  [17] plotly_4.10.4             httpuv_1.6.15
#  [19] sctransform_0.4.1         spam_2.11-1
#  [21] spatstat.sparse_3.1-0     reticulate_1.41.0.1
#  [23] cowplot_1.1.3             pbapply_1.7-2
#  [25] lubridate_1.9.4           zlibbioc_1.48.2
#  [27] abind_1.4-8               Rtsne_0.17
#  [29] purrr_1.0.4               R.utils_2.13.0
#  [31] BiocGenerics_0.48.1       tweenr_2.0.3
#  [33] circlize_0.4.16           IRanges_2.36.0
#  [35] S4Vectors_0.40.2          irlba_2.3.5.1
#  [37] listenv_0.9.1             spatstat.utils_3.1-2
#  [39] vegan_2.6-10              goftest_1.2-3
#  [41] RSpectra_0.16-2           dqrng_0.4.1
#  [43] spatstat.random_3.3-2     fitdistrplus_1.2-2
#  [45] parallelly_1.42.0         DelayedMatrixStats_1.24.0
#  [47] permute_0.9-7             codetools_0.2-20
#  [49] DelayedArray_0.28.0       ggforce_0.4.2
#  [51] tidyselect_1.2.1          shape_1.4.6.1
#  [53] farver_2.1.2              ScaledMatrix_1.10.0
#  [55] viridis_0.6.5             matrixStats_1.5.0
#  [57] stats4_4.3.3              spatstat.explore_3.3-4
#  [59] jsonlite_2.0.0            GetoptLong_1.0.5
#  [61] BiocNeighbors_1.20.2      tidygraph_1.3.1
#  [63] progressr_0.15.1          iterators_1.0.14
#  [65] ggridges_0.5.6            survival_3.8-3
#  [67] foreach_1.5.2             tools_4.3.3
#  [69] ica_1.0-3                 Rcpp_1.0.14
#  [71] glue_1.8.0                gridExtra_2.3
#  [73] SparseArray_1.2.4         mgcv_1.9-1
#  [75] MatrixGenerics_1.14.0     withr_3.0.2
#  [77] fastmap_1.2.0             rsvd_1.0.5
#  [79] digest_0.6.37             timechange_0.3.0
#  [81] R6_2.6.1                  mime_0.13
#  [83] ggprism_1.0.5             colorspace_2.1-1
#  [85] Cairo_1.6-2               scattermore_1.2
#  [87] tensor_1.5                spatstat.data_3.1-4
#  [89] R.methodsS3_1.8.2         generics_0.1.3
#  [91] data.table_1.17.0         graphlayouts_1.2.2
#  [93] httr_1.4.7                htmlwidgets_1.6.4
#  [95] S4Arrays_1.2.1            uwot_0.2.3
#  [97] pkgconfig_2.0.3           gtable_0.3.6
#  [99] ComplexHeatmap_2.18.0     lmtest_0.9-40
# [101] XVector_0.42.0            htmltools_0.5.8.1
# [103] dotCall64_1.2             clue_0.3-66
# [105] scales_1.4.0              png_0.1-8
# [107] spatstat.univar_3.1-2     snakecase_0.11.1
# [109] rjson_0.2.23              reshape2_1.4.4
# [111] checkmate_2.3.2           nlme_3.1-167
# [113] zoo_1.8-13                cachem_1.1.0
# [115] GlobalOptions_0.1.2       stringr_1.5.1
# [117] KernSmooth_2.23-26        parallel_4.3.3
# [119] miniUI_0.1.1.1            vipor_0.4.7
# [121] ggrastr_1.0.2             pillar_1.10.2
# [123] grid_4.3.3                vctrs_0.6.5
# [125] RANN_2.6.2                promises_1.3.2
# [127] BiocSingular_1.18.0       beachmat_2.18.1
# [129] xtable_1.8-4              cluster_2.1.8.1
# [131] beeswarm_0.4.0            paletteer_1.6.0
# [133] cli_3.6.5                 compiler_4.3.3
# [135] rlang_1.1.6               crayon_1.5.3
# [137] future.apply_1.11.3       labeling_0.4.3
# [139] rematch2_2.1.2            plyr_1.8.9
# [141] forcats_1.0.0             ggbeeswarm_0.7.2
# [143] stringi_1.8.7             viridisLite_0.4.2
# [145] deldir_2.0-4              BiocParallel_1.36.0
# [147] lazyeval_0.2.2            spatstat.geom_3.3-5
# [149] Matrix_1.6-5              RcppHNSW_0.6.0
# [151] patchwork_1.3.0           sparseMatrixStats_1.14.0
# [153] future_1.34.0             shiny_1.10.0
# [155] ROCR_1.0-11               igraph_2.1.4
# [157] memoise_2.0.1