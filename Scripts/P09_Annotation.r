.libPaths()
#[1] "/micromamba/envs/seurat+milo.new/lib/R/library"
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
library(findPC)
library(SingleCellExperiment)}


All3_integrated <- readRDS("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")

celltypes_subclustres02 <- c("5_0"="pDCs","5_1"="pDCs","4_2"= "NK cells","4_0"= "NK cells","4_1"= "NK cells","0_1"= "Monocytes","5_2"= "cDCs","0_4"= "Monocytes","0_3"= "Monocytes","0_0" = "Monocytes", "0_2"= "Monocytes","1_3"="Unknown", "1_1"= "B cells","1_2"= "B cells","1_0"= "B cells", "1_4"="ASCs","2_3"= "CD4+ ab T cells","2_2"= "CD4+ ab T cells", "2_0"= "CD4+ ab T cells","2_1"= "CD4+ ab T cells", "3_0"= "CD2- gd T cells","3_1"= "CD2- gd T cells","3_2"= "CD2- gd T cells","3_3"= "CD2- gd T cells")
Idents(All3_integrated)="sub.cluster.clustAll_res.0.2"
All3_integrated<- RenameIdents(All3_integrated, celltypes_subclustres02)
# Extract the current identities
celltypes_subclustres02 <- Idents(All3_integrated)
# Add the identities as a new column in the metadata
All3_integrated <- AddMetaData(All3_integrated, metadata = celltypes_subclustres02, col.name = "celltypes_subclustres02")
All3_integrated$celltypes_greek <- paste(All3_integrated$celltypes_subclustres02)  # Copy the "celltypes" column
All3_integrated$celltypes_greek <- gsub("gd", "γδ", All3_integrated$celltypes_greek)  # Replace "gd" with "γδ"

All3_integrated$celltypes_greek <- gsub("ab", "αβ", All3_integrated$celltypes_greek)  # Replace "ab" with "αβ"#table(colData(sce)$celltypes_greek)  # Create a frequency table
All3_integrated$DPI <- gsub("D(\\d+)", "\\1 DPI", All3_integrated$time_point)

saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")
saveRDS(All3_integrated,"/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered.rds")

All3_integrated2sce <- as.SingleCellExperiment(All3_integrated)
#add PCA & UMAP from Seurat object to SCE object
pca <- Embeddings(All3_integrated, reduction = "pca")
umap <- Embeddings(All3_integrated, reduction = "umap")
reducedDim(All3_integrated2sce, "PCA") <- pca
reducedDim(All3_integrated2sce, "UMAP") <- umap
saveRDS(All3_integrated2sce, file ="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered_sce.rds")

Jayne_NK_markers<- c("PAX5","CD19","CD79A","TCF4","SLA-DRB1","SLA-DRA","CR2","IRF8","CD86","TYROBP","HCST","KLRB1", "KLRK1","PRF1","CD6","CD5","CD8A","CD8B","CD4","TRDC","CD2","CD3E","CD3G","PRDM1","IRF4","CD93","CLEC12A","XBP1","FCER1A","FLT3","TLR4","NLRP3","CSF1R","CD163","CD14","SIRPA", "IL2RB" , "FCGR3A", "NCAM1", "CD69", "KLRD1", "COX6A2", "KIR2DL4", "NKG7", "NCR1", "ZMAT4")

custom_palette <- colorRampPalette(c("yellow", "red"))
DefaultAssay(All3_integrated)="RNA"
Idents(All3_integrated)="sub.cluster.clustAll_res.0.2"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclustres02_dot_plot_horizontal.pdf", width=12, height=8)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,x_lab_rotate=45,flip = TRUE,grid_color = "grey")
dev.off()

Idents(All3_integrated)="celltypes_subclustres02"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_celltypes_subclustres02_dot_plot_redyellowscale.pdf", width=9, height=10)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,grid_color = "grey",exp_color_min = -1,exp_color_max = 2, colors_use_exp = custom_palette(100),x_lab_rotate=90)
dev.off()

Idents(All3_integrated)="celltypes_subclustres02"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_celltypes_dot_plot.pdf", width=9, height=10)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,grid_color = "grey")
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_celltypes_dot_plot_horizontal.pdf", width=12, height=4)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,x_lab_rotate=45,flip = TRUE,grid_color = "grey")
dev.off()



celltype_palette <- brewer.pal(n = 12, name = "Set3")
# I know #ffffb3 is the super light yellow, so replace it with "#FFED6F"
celltype_palette2 <- celltype_palette
celltype_palette2[2] <- "#FFED6F"

library(Cairo)
Idents(All3_integrated)="celltypes_greek"
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Greekcelltypes_dot_plot_horizontal.pdf", width=12, height=4)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,colors_use_idents=celltype_palette2,x_lab_rotate=45,flip = TRUE,grid_color = "grey")
dev.off()

All3_integrated$cluster_celltype <- as.factor(paste0((All3_integrated$sub.cluster.clustAll_res.0.2), " ", All3_integrated$celltypes_greek))
Idents(All3_integrated)="cluster_celltype"
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_ClusterGreekcelltypes_dot_plot_horizontal.pdf", width=12, height=8)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,cluster_feature = FALSE,cluster_ident = FALSE,x_lab_rotate=45,flip = TRUE)
dev.off()
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Clustercelltypes_dot_plot_horizontal.pdf", width=12, height=8)
Clustered_DotPlot(All3_integrated, features=Jayne_NK_markers,cluster_feature = FALSE,cluster_ident = FALSE,x_lab_rotate=45,flip = TRUE)
dev.off()



Idents(All3_integrated)="celltypes_greek"
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_celltypes_UMAP.pdf", width=10, height=10)
DimPlot_scCustom(All3_integrated,colors_use = celltype_palette2) +theme(aspect.ratio = 1) + ggtitle("STACAs integration via Animal, ScranNorm_Regressed_PCA_d10_1284HVF")
dev.off()

#remove "celltypes_greek"== "Unknown" and plot again
All3_integrated_noUnknown <- subset(All3_integrated, subset = celltypes_greek != "Unknown")
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Greekcelltypes_dot_plot_horizontal_paper.pdf", width=12, height=4)
Clustered_DotPlot(All3_integrated_noUnknown, features=Jayne_NK_markers,colors_use_idents=celltype_palette2,x_lab_rotate=45,flip = TRUE,show_parent_dend_line = FALSE)
dev.off()

CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_celltypes_UMAP_paper.pdf", width=10, height=10)
DimPlot_scCustom(All3_integrated_noUnknown,colors_use = celltype_palette2,label = FALSE,raster = TRUE) +theme(aspect.ratio = 1) 
dev.off()

Idents(All3_integrated)<-"sub.cluster.clustAll_res.0.2"
cellproportions2<-prop.table(table(Idents(All3_integrated), All3_integrated$time_point), margin = 2)
cellproportionsdfs2<- write.table(cellproportions2, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_time_prop.txt", sep="\t", quote=F)
cluster_prop_time<-as.data.frame(cellproportions2)
colnames(cluster_prop_time) <-c("cluster","DPI", "Proportion")
cluster_prop_time$cluster<-as.factor(cluster_prop_time$cluster)
#set the order of cluster
cluster_prop_time$cluster <- factor(cluster_prop_time$cluster, levels=c("0_0", "0_1", "0_2","0_3","0_4", "1_0", "1_1","1_2","1_3","1_4", "2_0", "2_1", "2_2", "2_3","3_0", "3_1","3_2","3_3", "4_0", "4_1","4_2", "5_0", "5_1","5_2"))

#Make bar plot fill by day cluster = 100%
 # Color-blind friendly colors
colors <- c("#D55E00", "#0072B2", "#009E73") 
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Barplot_time_cluster.pdf",width=20, height=10)
plot1<-ggplot(cluster_prop_time, aes(x = cluster, y = Proportion, fill = DPI)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) + # Use position dodge for unstacked bars
  theme_minimal() +
  labs(title = "Proportions of cells across the whole dataset",
       x = "Cluster",
       y = "Proportion") +
  theme(axis.title = element_blank(), # Remove axis titles,
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,color = "black",size = 25,face = "bold"),
        axis.text.y = element_text(color = "black",size = 25,face = "bold"),
        legend.title = element_blank(), # Remove legend title
        legend.text = element_text(size = 25, face = "bold"),
        legend.key.size = unit(1.5, "cm"),  # Set size of legend squares
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors) # Use the same color palette
plot1
dev.off()


Idents(All3_integrated)<-"DPI"
cellproportions4<-prop.table(table(Idents(All3_integrated),All3_integrated$sub.cluster.clustAll_res.0.2), margin = 2)
cellproportionsdfs4<- write.table(cellproportions4, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_time_by_cluster_prop.txt", sep="\t", quote=F)
#make rows as time points and columns as clusters
cellproportions4 <- as.data.frame(cellproportions4) %>%
  as_tibble()
colnames(cellproportions4) <-c("DPI","SubCluster", "Proportion")
cellproportions4$SubCluster<-as.factor(cellproportions4$SubCluster)
#set the order of cluster
cellproportions4$SubCluster <- factor(cellproportions4$SubCluster, levels=c("0_0", "0_1", "0_2","0_3","0_4", "1_0", "1_1","1_2","1_3","1_4", "2_0", "2_1", "2_2", "2_3","3_0", "3_1","3_2","3_3", "4_0", "4_1","4_2", "5_0", "5_1","5_2"))
#Make bar plot fill by day cluster = 100%
 # Color-blind friendly colors
colors <- c("#D55E00", "#0072B2", "#009E73") 
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Barplot_cluster_by_time.pdf",width=20, height=10)
plot1<-ggplot(cellproportions4, aes(x = SubCluster, y = Proportion, fill = DPI)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) + # Use position dodge for unstacked bars
  theme_minimal() +
  labs(title = "Proportions of cells within SubCluster",
       x = "SubCluster",
       y = "Proportion") +
  theme(axis.title = element_blank(), # Remove axis titles,
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,color = "black",size = 25,face = "bold"),
        axis.text.y = element_text(color = "black",size = 25,face = "bold"),
        legend.title = element_blank(), # Remove legend title
        legend.text = element_text(size = 25, face = "bold"),
        legend.key.size = unit(1.5, "cm"),  # Set size of legend squares
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors) # Use the same color palette
plot1
dev.off()


Idents(All3_integrated)<-"Animal"
cellproportions3<-prop.table(table(Idents(All3_integrated), All3_integrated$sub.cluster.clustAll_res.0.2), margin = 2)
cellproportionsdfs3<- write.table(cellproportions3, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_Animal.txt", sep="\t", quote=F)
cluster_prop_Animal<-as.data.frame(cellproportions3)
colnames(cluster_prop_Animal) <-c("Animal","SubCluster", "Proportion")
cluster_prop_Animal$SubCluster<-as.factor(cluster_prop_Animal$SubCluster)
cluster_prop_Animal$SubCluster <- factor(cluster_prop_Animal$SubCluster, levels=c("0_0", "0_1", "0_2","0_3","0_4", "1_0", "1_1","1_2","1_3","1_4", "2_0", "2_1", "2_2", "2_3","3_0", "3_1","3_2","3_3", "4_0", "4_1","4_2", "5_0", "5_1","5_2"))

#Make starcked bar plot fill by day cluster = 100%
 # Color-blind friendly colors
colors2 <- c("#E69F00",  "#F0E442",  "#CC79A7", "#009E73", "#0072B2")
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Barplot_time_Animal.pdf",width=20, height=10)
plot1<-ggplot(cluster_prop_Animal,  aes(x = SubCluster, y = Proportion, fill = Animal)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) + # Use position dodge for unstacked bars
  theme_minimal() +
  labs(title = "Proportions of cells from animal within SubCluster",
       x = "SubCluster",
       y = "Proportion") +
  theme(axis.title = element_blank(), # Remove axis titles,
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,color = "black",size = 25,face = "bold"),
        axis.text.y = element_text(color = "black",size = 25,face = "bold"),
        legend.title = element_blank(), # Remove legend title
        legend.text = element_text(size = 25, face = "bold"),
        legend.key.size = unit(1.5, "cm"),  # Set size of legend squares
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors2)
plot1  
dev.off()

# number of cells per celltype
cellnumdf<-as.data.frame(table(All3_integrated$celltypes_subclustres02))
colnames(cellnumdf)<-c("celltype","cells")
cellnumdfs<- write.csv(cellnumdf, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_Ncelltypes.csv")
Idents(All3_integrated)="celltypes_subclustres02"
celltype_time<-table(Idents(All3_integrated), All3_integrated$time_point)
write.table(celltype_time, file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_time_Number.txt", sep="\t", quote=F)
celltype_time_prop<-prop.table(table(Idents(All3_integrated), All3_integrated$time_point), margin = 2)
write.table(celltype_time_prop, file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_time_Prop.txt", sep="\t", quote=F)
All3_integrated$Animal_time<-paste0(All3_integrated$Animal, "_", All3_integrated$time_point)
celltype_Animal_time<-table(Idents(All3_integrated), All3_integrated$Animal_time)
write.table(celltype_Animal_time, file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Animal_Number_Prop.txt", sep="\t", quote=F)

Idents(All3_integrated)="Animal_time"
celltype_Animal_time_prop<-prop.table(table(Idents(All3_integrated), All3_integrated$celltypes_subclustres02), margin = 2)
write.table(celltype_Animal_time_prop, file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Animal_time_Prop.txt", sep="\t", quote=F)

# Get the colorblind_palett for DPI
colorblind_palette <- ColorBlind_Pal()
# Get the dark2 palette for Animal
dark2_palette <- Dark2_Pal()
#Make polychrome colorpalette for AAI Celltypes
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 12, palette = "polychrome", shuffle_pal = FALSE)
# Remove the second color from the palette
polychrome_pal <- polychrome_pal[-2]

CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_4_UMAPs.pdf", width=40, height=10)
P1<-DimPlot_scCustom(All3_integrated, group.by = "celltypes_greek",colors_use=celltype_palette2,repel=TRUE,raster = TRUE)
P2<-DimPlot_scCustom(All3_integrated, group.by = "DPI",colors_use = colorblind_palette,raster = TRUE)
P3<-DimPlot_scCustom(All3_integrated, group.by = "sub.cluster.clustAll_res.0.2",colors_use = DiscretePalette_scCustomize(num_colors = 43,palette = "varibow"),raster = TRUE)
P4<-DimPlot_scCustom(All3_integrated, group.by = "Animal",colors_use = dark2_palette,raster = TRUE)
P3|P1|P2|P4
dev.off()
# Make UMAP colored by Day & split by day
colorblind_palette <- ColorBlind_Pal()
pdf(file="//scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_annotated_split_Time.pdf", width=30, height=10)
Idents(All3_integrated)<-"DPI"
DimPlot_scCustom(All3_integrated,split.by="DPI",colors_use = colorblind_palette,label=FALSE,raster =TRUE)
dev.off()
# Make UMAP colored by Animal & split by Animal
colorblind_palette <- ColorBlind_Pal()
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_annotated_split_Animal.pdf", width=30, height=10)
Idents(All3_integrated)<-"Animal"
DimPlot_scCustom(All3_integrated,split.by="Animal",colors_use = colorblind_palette,label=FALSE,raster =TRUE)
dev.off()

colorblind_palette <- ColorBlind_Pal()
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_annotated_celltype_split_time_UMAP.pdf", width=30, height=10)
Idents(All3_integrated)<-"celltypes_greek"
DimPlot_scCustom(All3_integrated,split.by="DPI",colors_use = celltype_palette2,label=FALSE,raster =TRUE)
dev.off()



Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.05"
NKcells <- subset(All3_integrated, idents=c("4_0"))
Idents(NKcells) <- "sub.cluster.clustAll_res.0.09"
#remove reslutions above 0.125
{NKcells$"sub.cluster.clust4_res.0.15" <- NULL
NKcells$"sub.cluster.clust4_res.0.1525" <- NULL
NKcells$"sub.cluster.clust4_res.0.155" <- NULL
NKcells$"sub.cluster.clust4_res.0.16" <- NULL
NKcells$"sub.cluster.clust4_res.0.175" <- NULL
NKcells$"sub.cluster.clust4_res.0.2" <- NULL}
pdf("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclustered_NKcells_clustree.pdf", width=6, height=8)
clustree(NKcells,  prefix = "sub.cluster.clust4_res.")
dev.off()

# Replace the "sub.cluster.clustAll_res.0.09"==4_0 with "NKC0", and sub.cluster.clustAll_res.0.09"==4_1 with "NKC1"
NKcells$subcluster_name <-  NKcells$sub.cluster.clustAll_res.0.09
NKcells$subcluster_name <- gsub("4_0", "NKC0", NKcells$subcluster_name)
NKcells$subcluster_name <- gsub("4_1", "NKC1", NKcells$subcluster_name)

Idents(NKcells) <- "subcluster_name"
DefaultAssay(NKcells) <- "RNA"
NKmarkers <- c("FCGR3A","FCER1G","ITGAM","PRF1","GZMB","GZMK","GNLY","IL7R","NKG7","CCL4","XCL2","NCAM1","SELL","B3GAT1","CD2","LAG3","KLRC1","KLRC2","CD27","MKI67","TOP2A","E2F1")
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_NKcellsOnly_dot_plot.pdf")
DotPlot_scCustom(NKcells, features=NKmarkers,x_lab_rotate = TRUE,)
dev.off()

cellproportions2<-prop.table(table( NKcells$DPI, NKcells$subcluster_name), margin = 2)
cluster_prop_time<-as.data.frame(cellproportions2)
colnames(cluster_prop_time) <-c("DPI","subcluster", "Proportion")
cluster_prop_time$subcluster<-as.factor(cluster_prop_time$subcluster)
head(cluster_prop_time)

#Make starcked bar plot fill by day cluster = 100%
 # Color-blind friendly colors
colors <- c("#D55E00", "#0072B2", "#009E73") 


cellnumbers2<-table(NKcells$DPI, NKcells$subcluster_name)
write.table(cellnumbers2, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/NKcells_FindSubCluster_Subcluster_res009_Cellnumbers_time_cluster100.txt", sep="\t", quote=F)

pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/NKcells_FindSubCluster_Subcluster_res009_Barplot_time_cluster100.pdf")
plot2 <- ggplot(cluster_prop_time, aes(x = subcluster, y = Proportion, fill = DPI)) +
  geom_bar(stat = "identity", position = "dodge") + # Use position dodge for unstacked bars
  theme_minimal() +
  labs(title = "Proportions of cells from time points across NKcells subclusters res0.09, cluster = 100%",
       x = NULL, #"Subcluster",
       y = "Proportion of Cluster") +
theme(
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(margin = margin(r = 15), size = 18, face = "bold", color = "black"), # <-- add space here (15 pt)
    axis.text = element_text(color = "black", size = 18, face = "bold"),
    axis.text.x = element_text( hjust = 0.5, vjust = 0.5),#angle = 315,
    legend.title =  element_blank(),
    legend.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.position = "top"
  ) +  scale_fill_manual(values = colors) # Use the same color palette
plot2
dev.off()

Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.05"
Monocytes <- subset(All3_integrated, idents=c("0_0"))
Idents(Monocytes) <- "sub.cluster.clustAll_res.0.09"
#remove reslutions above 0.125
{Monocytes$"sub.cluster.clust0_res.0.15" <- NULL
Monocytes$"sub.cluster.clust0_res.0.1525" <- NULL
Monocytes$"sub.cluster.clust0_res.0.155" <- NULL
Monocytes$"sub.cluster.clust0_res.0.16" <- NULL
Monocytes$"sub.cluster.clust0_res.0.175" <- NULL
Monocytes$"sub.cluster.clust0_res.0.2" <- NULL}
pdf("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_subclustered_Monocytes_clustree.pdf", width=6, height=8)
clustree(Monocytes,  prefix = "sub.cluster.clust0_res.")
dev.off()

# Replace the "sub.cluster.clustAll_res.0.09"==0_0 with "MC0", and sub.cluster.clustAll_res.0.09"==0_1 with "MC1"
Monocytes$subcluster_name <- Monocytes$sub.cluster.clustAll_res.0.09
Monocytes$subcluster_name <- gsub("0_0", "MC0", Monocytes$subcluster_name)
Monocytes$subcluster_name <- gsub("0_1", "MC1", Monocytes$subcluster_name)

cellproportions2<-prop.table(table( Monocytes$DPI, Monocytes$subcluster_name), margin = 2)
cluster_prop_time<-as.data.frame(cellproportions2)
colnames(cluster_prop_time) <-c("DPI","subcluster", "Proportion")
cluster_prop_time$subcluster<-as.factor(cluster_prop_time$subcluster)
head(cluster_prop_time)

#Make starcked bar plot fill by day cluster = 100%
 # Color-blind friendly colors
colors <- c("#D55E00", "#0072B2", "#009E73") 


cellnumbers2<-table(Monocytes$DPI, Monocytes$subcluster_name)
write.table(cellnumbers2, file = "/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/Monocytes_FindSubCluster_Subcluster_res009_Cellnumbers_time_cluster100.txt", sep="\t", quote=F)

pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/Monocytes_FindSubCluster_Subcluster_res009_Barplot_time_cluster100.pdf")
plot2 <- ggplot(cluster_prop_time, aes(x = subcluster, y = Proportion, fill = DPI)) +
  geom_bar(stat = "identity", position = "dodge") + # Use position dodge for unstacked bars
  theme_minimal() +
  labs(title = "Proportions of cells from time points across Monocytes subclusters res0.09, cluster = 100%",
       x = NULL, #"Subcluster",
       y = "Proportion of Cluster") +
theme(
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(margin = margin(r = 15), size = 18, face = "bold", color = "black"), # <-- add space here (15 pt)
    axis.text = element_text(color = "black", size = 18, face = "bold"),
    axis.text.x = element_text( hjust = 0.5, vjust = 0.5), #angle = 315,
    legend.title =  element_blank(),
    legend.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0, size = 10, face = "bold"),
    legend.position = "top"
  ) +  scale_fill_manual(values = colors) # Use the same color palette
plot2
dev.off()

colorblind_palette <- ColorBlind_Pal()
Idents(All3_integrated) <- "sub.cluster.clustAll_res.0.09"
pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_highlight_Monocyte_clusters.pdf", width=10, height=10)
Cluster_Highlight_Plot(seurat_object = All3_integrated, cluster_name = c("0_0","0_1"), highlight_color = c(colorblind_palette),raster=TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_highlight_NKcells_clusters.pdf", width=10, height=10)
Cluster_Highlight_Plot(seurat_object = All3_integrated, cluster_name = c("4_0","4_1"), highlight_color = c(colorblind_palette),raster=TRUE)
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
#   [1] RcppAnnoy_0.0.22            splines_4.3.3
#   [3] later_1.4.2                 prismatic_1.1.2
#   [5] bitops_1.0-9                tibble_3.2.1
#   [7] R.oo_1.27.0                 polyclip_1.10-7
#   [9] janitor_2.2.1               fastDummies_1.7.5
#  [11] lifecycle_1.0.4             doParallel_1.0.17
#  [13] globals_0.16.3              lattice_0.22-6
#  [15] MASS_7.3-60                 magrittr_2.0.3
#  [17] plotly_4.10.4               httpuv_1.6.15
#  [19] sctransform_0.4.1           spam_2.11-1
#  [21] spatstat.sparse_3.1-0       reticulate_1.41.0.1
#  [23] cowplot_1.1.3               pbapply_1.7-2
#  [25] lubridate_1.9.4             zlibbioc_1.48.2
#  [27] abind_1.4-8                 GenomicRanges_1.54.1
#  [29] Rtsne_0.17                  purrr_1.0.4
#  [31] R.utils_2.13.0              RCurl_1.98-1.17
#  [33] BiocGenerics_0.48.1         tweenr_2.0.3
#  [35] GenomeInfoDbData_1.2.11     circlize_0.4.16
#  [37] IRanges_2.36.0              S4Vectors_0.40.2
#  [39] irlba_2.3.5.1               listenv_0.9.1
#  [41] spatstat.utils_3.1-2        vegan_2.6-10
#  [43] goftest_1.2-3               RSpectra_0.16-2
#  [45] dqrng_0.4.1                 spatstat.random_3.3-2
#  [47] fitdistrplus_1.2-2          parallelly_1.42.0
#  [49] DelayedMatrixStats_1.24.0   permute_0.9-7
#  [51] codetools_0.2-20            DelayedArray_0.28.0
#  [53] ggforce_0.4.2               tidyselect_1.2.1
#  [55] shape_1.4.6.1               farver_2.1.2
#  [57] ScaledMatrix_1.10.0         viridis_0.6.5
#  [59] matrixStats_1.5.0           stats4_4.3.3
#  [61] spatstat.explore_3.3-4      jsonlite_2.0.0
#  [63] GetoptLong_1.0.5            BiocNeighbors_1.20.2
#  [65] tidygraph_1.3.1             progressr_0.15.1
#  [67] iterators_1.0.14            ggridges_0.5.6
#  [69] survival_3.8-3              foreach_1.5.2
#  [71] tools_4.3.3                 ica_1.0-3
#  [73] Rcpp_1.0.14                 glue_1.8.0
#  [75] gridExtra_2.3               SparseArray_1.2.4
#  [77] mgcv_1.9-1                  MatrixGenerics_1.14.0
#  [79] GenomeInfoDb_1.38.8         withr_3.0.2
#  [81] fastmap_1.2.0               rsvd_1.0.5
#  [83] digest_0.6.37               timechange_0.3.0
#  [85] R6_2.6.1                    mime_0.13
#  [87] ggprism_1.0.5               colorspace_2.1-1
#  [89] Cairo_1.6-2                 scattermore_1.2
#  [91] tensor_1.5                  spatstat.data_3.1-4
#  [93] R.methodsS3_1.8.2           utf8_1.2.4
#  [95] generics_0.1.3              data.table_1.17.0
#  [97] graphlayouts_1.2.2          httr_1.4.7
#  [99] htmlwidgets_1.6.4           S4Arrays_1.2.1
# [101] uwot_0.2.3                  pkgconfig_2.0.3
# [103] gtable_0.3.6                ComplexHeatmap_2.18.0
# [105] lmtest_0.9-40               SingleCellExperiment_1.24.0
# [107] XVector_0.42.0              htmltools_0.5.8.1
# [109] dotCall64_1.2               clue_0.3-66
# [111] Biobase_2.62.0              scales_1.4.0
# [113] png_0.1-8                   spatstat.univar_3.1-2
# [115] snakecase_0.11.1            rjson_0.2.23
# [117] reshape2_1.4.4              nlme_3.1-167
# [119] zoo_1.8-13                  cachem_1.1.0
# [121] GlobalOptions_0.1.2         stringr_1.5.1
# [123] KernSmooth_2.23-26          parallel_4.3.3
# [125] miniUI_0.1.1.1              vipor_0.4.7
# [127] ggrastr_1.0.2               pillar_1.10.2
# [129] grid_4.3.3                  vctrs_0.6.5
# [131] RANN_2.6.2                  promises_1.3.2
# [133] BiocSingular_1.18.0         beachmat_2.18.1
# [135] xtable_1.8-4                cluster_2.1.8.1
# [137] beeswarm_0.4.0              paletteer_1.6.0
# [139] cli_3.6.5                   compiler_4.3.3
# [141] rlang_1.1.6                 crayon_1.5.3
# [143] future.apply_1.11.3         labeling_0.4.3
# [145] rematch2_2.1.2              plyr_1.8.9
# [147] forcats_1.0.0               ggbeeswarm_0.7.2
# [149] stringi_1.8.7               viridisLite_0.4.2
# [151] deldir_2.0-4                BiocParallel_1.36.0
# [153] lazyeval_0.2.2              spatstat.geom_3.3-5
# [155] Matrix_1.6-5                RcppHNSW_0.6.0
# [157] patchwork_1.3.0             sparseMatrixStats_1.14.0
# [159] future_1.34.0               shiny_1.10.0
# [161] SummarizedExperiment_1.32.0 ROCR_1.0-11
# [163] igraph_2.1.4                memoise_2.0.1
