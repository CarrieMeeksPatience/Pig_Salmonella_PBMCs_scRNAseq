.libPaths()
#[1] "/micromamba/envs/miloR2.0/lib/R/library"
{set.seed(123)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scRNAseq)
library(scuttle)
library(irlba)
library(BiocParallel)
library(ggplot2)
library(cowplot)
library(Cairo)}



sce <- readRDS("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered_sce.rds")
#reducedDims(sce)
#Add animal_time to metadata
colData(sce)$Animal_time <- paste(colData(sce)$Animal, colData(sce)$time_point, sep="_")
colData(sce)$Animal_time<-as.character(colData(sce)$Animal_time)
saveRDS(sce,file="/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered_sce.rds")

sce <- readRDS("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_2026_04_16_subclustered_sce.rds")
head(colData(sce))
#remove cells where celltypes_greek=="Unknown"
sce <- sce[, sce$celltypes_greek != "Unknown"]
#Split by time points
sce_D20 <- sce[, sce$time_point %in% c("D2", "D0")]

#Make milo object from the SCE object
D20milo <- Milo(sce_D20) 


#Add UMAP from SCE object to milo object
reducedDim(D20milo, "UMAP") <- reducedDim(sce_D20, "UMAP")


####=======  k=# of nearest neighbors to include(# of samples x 5 = 75) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 5=75), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D20milo_k75 <- buildGraph(D20milo, k = 75, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D20milo_k75 <- makeNhoods(D20milo_k75, prop = 0.1, k = 75, d = 10, refined = TRUE, reduced_dims = "UMAP",refinement_scheme="graph") 


#add animal_time to metadata
colData(D20milo_k75)$ObsID <- paste(colData(D20milo_k75)$animal, colData(D20milo_k75)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D20milo_k75 <- countCells(D20milo_k75, meta.data = data.frame(colData(D20milo_k75)), samples="Animal_time") 
D20milo_k75

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k75Ncounts<-nhoodCounts(D20milo_k75)
k75Ncounts_s<-saveRDS(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D2v0_2026_04_17.rds")
k75Ncountsm<-as.matrix(k75Ncounts)
k75Ncountsdf<-as.data.frame(k75Ncountsm)
k75Ncounts<-k75Ncountsdf
k75Ncounts_s<-write.csv(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D2v0_2026_04_17.csv")
avg_k75Ncounts<-colMeans(k75Ncounts)
print(avg_k75Ncounts)
median_k75Ncounts<-sapply(k75Ncounts,median)
print(median_k75Ncounts)
IQR_k75Ncounts<-sapply(k75Ncounts,IQR)
print(IQR_k75Ncounts)
save<-saveRDS(D20milo_k75,file="/scRNAseq/Sal_5pigs_2026/DA/milo_D2v0_k75_pca10_d10_prop0.1_2026_04_17.rds")
####=======  k=# of nearest neighbors to include(# of samples x 3 = 45) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 3=45), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D20milo_k45 <- buildGraph(D20milo, k = 45, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D20milo_k45 <- makeNhoods(D20milo_k45, prop = 0.1, k = 45, d = 10, refined = TRUE, reduced_dims = "UMAP",refinement_scheme="graph") 


#add animal_time to metadata
colData(D20milo_k45)$ObsID <- paste(colData(D20milo_k45)$animal, colData(D20milo_k45)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D20milo_k45 <- countCells(D20milo_k45, meta.data = data.frame(colData(D20milo_k45)), samples="Animal_time") 
D20milo_k45

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k45Ncounts<-nhoodCounts(D20milo_k45)
k45Ncounts_s<-saveRDS(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/K45Ncounts_D2v0_2026_04_17.rds")
k45Ncountsm<-as.matrix(k45Ncounts)
k45Ncountsdf<-as.data.frame(k45Ncountsm)
k45Ncounts<-k45Ncountsdf
k45Ncounts_s<-write.csv(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/K45Ncounts_D2v0_2026_04_17.csv")
avg_k45Ncounts<-colMeans(k45Ncounts)
print(avg_k45Ncounts)
median_k45Ncounts<-sapply(k45Ncounts,median)
print(median_k45Ncounts)
IQR_k45Ncounts<-sapply(k45Ncounts,IQR)
print(IQR_k45Ncounts)

det_K<-rbind(avg_k75Ncounts,median_k75Ncounts,IQR_k75Ncounts,avg_k45Ncounts,median_k45Ncounts,IQR_k45Ncounts)

det_k_s<-write.csv(det_K, file="/scRNAseq/Sal_5pigs_2026/DA/k45/K45_75_D2v0_avg_IQR_2026_04_17.csv")

k45Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D20milo_k45)),
                        "Var"=rowVars(nhoodCounts(D20milo_k45)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_k45_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k45Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

k75Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D20milo_k75)),
                        "Var"=rowVars(nhoodCounts(D20milo_k75)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k75Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

#note by Carrie: counts= # of neighbourhoods neighbourhood size = # of cells in neighbourhood
pdf("/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D20milo_k75)
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_k45_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D20milo_k45)
dev.off()


#Design model matrix from metadata, make sure the rownames match with nhoodCounts column names & are in same order
pbmc.design <- data.frame(colData(D20milo_k75))[,c("Animal_time", "time_point", "Animal")]
pbmc.design$Animal <- as.factor(pbmc.design$Animal) 
pbmc.design <- distinct(pbmc.design)
rownames(pbmc.design) <- pbmc.design$Animal_time
pbmc.design <- pbmc.design[colnames(nhoodCounts(D20milo_k75)), , drop=FALSE]
head(pbmc.design)

#Calculate Nhood distance to later use to account for the overlap between neighbourhoods & build Nhood graph
D20milo_k75 <- calcNhoodDistance(D20milo_k75, d = 10)
#could try D20milo_k75 <- buildNhoodGraph(D20milo_k75, overlap = 20)
D20milo_k75 <- buildNhoodGraph(D20milo_k75, overlap = 60)
save<-saveRDS(D20milo_k75, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_pca10_d10_prop0.1_2026_04_17.rds")


D20milo_k45 <- calcNhoodDistance(D20milo_k45, d = 10)
D20milo_k45 <- buildNhoodGraph(D20milo_k45, overlap = 30)
save<-saveRDS(D20milo_k45, file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_k45_pca10_d10_prop0.1_2026_04_17.rds")


#Run DA glmm
##run tasks in a serial manner
bpparam <- SerialParam()
register(bpparam)
##use non-negative least squares (NNLS) Haseman-Elston solver,REML is Restricted Maximum Likelihood),fdr.weighting=spatial FDR weighting scheme to use
print("D20 k75")
glm_results_D20 <- testNhoods(D20milo_k75, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)
which(checkSeparation(D20milo_k75, design.df=pbmc.design, condition="time_point", min.val=10))


da_results_D20 <- testNhoods(D20milo_k75, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)



comp.da <- merge(da_results_D20, glm_results_D20, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()



## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_NhoodGraph_2026_04_17.pdf", width=10, height=10)
plotNhoodGraphDA(D20milo_k75, da_results_D20, layout="UMAP",alpha=0.1) 
dev.off()
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_NhoodGraph_2026_04_17_v2.pdf", width=10, height=10)
plotNhoodGraphDA(D20milo_k75, da_results_D20, layout="UMAP",alpha=0.1,node_stroke=0.1) 
dev.off()

# Add cell types to results 
da_results_D20 <- annotateNhoods(D20milo_k75, da_results_D20, coldata_col = "celltypes_greek")
head(da_results_D20)

saveRDS(da_results_D20, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_2026_04_17.rds")
write.table(da_results_D20, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_2026_04_17.txt", sep="\t", row.names=FALSE)

rm(comp.da)
### Try with K45
print("D20 k45")

which(checkSeparation(D20milo_k45, design.df=pbmc.design, condition="time_point", min.val=10))

glm_results_D20_k45 <- testNhoods(D20milo_k45, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)

da_results_D20_k45 <- testNhoods(D20milo_k45, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)

comp.da <- merge(da_results_D20_k45, glm_results_D20_k45, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_k45_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()

## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_45_DA_glmm_NhoodGraph_2026_04_17.pdf")
plotNhoodGraphDA(D20milo_k45,da_results_D20_k45, layout="UMAP",alpha=0.1) 
dev.off()

# Add cell types to results 
da_results_D20_k45 <- annotateNhoods(D20milo_k45, da_results_D20_k45, coldata_col = "celltypes_greek")
head(da_results_D20_k45)

saveRDS(da_results_D20_k45, file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D20_k45_DA_glmm_2026_04_17.rds")


print("Run D8 vs D0")
#Split by time points
sce_D80 <- sce[, sce$time_point %in% c("D0", "D8")]


#Make milo object from the SCE object
D80milo <- Milo(sce_D80) 


#Add UMAP from SCE object to milo object
reducedDim(D80milo, "UMAP") <- reducedDim(sce_D80, "UMAP")


####=======  k=# of nearest neighbors to include(# of samples x 5 = 75) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 5=75), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D80milo_k75 <- buildGraph(D80milo, k = 75, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D80milo_k75 <- makeNhoods(D80milo_k75, prop = 0.1, k = 75, d = 10, refined = TRUE, reduced_dims = "UMAP",refinement_scheme="graph") 

#add animal_time to metadata
colData(D80milo_k75)$ObsID <- paste(colData(D80milo_k75)$animal, colData(D80milo_k75)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D80milo_k75 <- countCells(D80milo_k75, meta.data = data.frame(colData(D80milo_k75)), samples="Animal_time") 
D80milo_k75

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k75Ncounts<-nhoodCounts(D80milo_k75)
k75Ncounts_s<-saveRDS(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D8v0_2026_04_17.rds")
k75Ncountsm<-as.matrix(k75Ncounts)
k75Ncountsdf<-as.data.frame(k75Ncountsm)
k75Ncounts<-k75Ncountsdf
k75Ncounts_s<-write.csv(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D8v0_2026_04_17.csv")
avg_k75Ncounts<-colMeans(k75Ncounts)
print(avg_k75Ncounts)
median_k75Ncounts<-sapply(k75Ncounts,median)
print(median_k75Ncounts)
IQR_k75Ncounts<-sapply(k75Ncounts,IQR)
print(IQR_k75Ncounts)
save<-saveRDS(D80milo_k75,file="/scRNAseq/Sal_5pigs_2026/DA/milo_D8v0_k75_pca10_d10_prop0.1_2026_04_17.rds")
####=======  k=# of nearest neighbors to include(# of samples x 3 = 45) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 3=45), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D80milo_k45 <- buildGraph(D80milo, k = 45, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D80milo_k45 <- makeNhoods(D80milo_k45, prop = 0.1, k = 45, d = 10, refined = TRUE,reduced_dims = "UMAP", refinement_scheme="graph") 

#add animal_time to metadata
colData(D80milo_k45)$ObsID <- paste(colData(D80milo_k45)$animal, colData(D80milo_k45)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D80milo_k45 <- countCells(D80milo_k45, meta.data = data.frame(colData(D80milo_k45)), samples="Animal_time") 
D80milo_k45

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k45Ncounts<-nhoodCounts(D80milo_k45)
k45Ncounts_s<-saveRDS(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/k45Ncounts_D8v0_2026_04_17.rds")
k45Ncountsm<-as.matrix(k45Ncounts)
k45Ncountsdf<-as.data.frame(k45Ncountsm)
k45Ncounts<-k45Ncountsdf
k45Ncounts_s<-write.csv(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/k45Ncounts_D8v0_2026_04_17.csv")
avg_k45Ncounts<-colMeans(k45Ncounts)
print(avg_k45Ncounts)
median_k45Ncounts<-sapply(k45Ncounts,median)
print(median_k45Ncounts)
IQR_k45Ncounts<-sapply(k45Ncounts,IQR)
print(IQR_k45Ncounts)

det_K<-rbind(avg_k75Ncounts,median_k75Ncounts,IQR_k75Ncounts,avg_k45Ncounts,median_k45Ncounts,IQR_k45Ncounts)

det_k_s<-write.csv(det_K, file="/scRNAseq/Sal_5pigs_2026/DA/K45_75_D8v0_avg_IQR_2026_04_17.csv")

k45Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D80milo_k45)),
                        "Var"=rowVars(nhoodCounts(D80milo_k45)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D80_k45_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k45Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

k75Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D80milo_k75)),
                        "Var"=rowVars(nhoodCounts(D80milo_k75)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k75Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

#note by Carrie: counts= # of neighbourhoods neighbourhood size = # of cells in neighbourhood
pdf("/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D80milo_k75)
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D80_k45_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D80milo_k45)
dev.off()


#Design model matrix from metadata, make sure the rownames match with nhoodCounts column names & are in same order
pbmc.design <- data.frame(colData(D80milo_k75))[,c("Animal_time", "time_point", "Animal")]
pbmc.design$Animal <- as.factor(pbmc.design$Animal) 
pbmc.design <- distinct(pbmc.design)
rownames(pbmc.design) <- pbmc.design$Animal_time
pbmc.design <- pbmc.design[colnames(nhoodCounts(D80milo_k75)), , drop=FALSE]
head(pbmc.design)

#Calculate Nhood distance to later use to account for the overlap between neighbourhoods & build Nhood graph
D80milo_k75 <- calcNhoodDistance(D80milo_k75, d = 10)
D80milo_k75 <- buildNhoodGraph(D80milo_k75,overlap=60)
save<-saveRDS(D80milo_k75, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_pca10_d10_prop0.1_2026_04_17.rds")


D80milo_k45 <- calcNhoodDistance(D80milo_k45, d = 10)
D80milo_k45 <- buildNhoodGraph(D80milo_k45,overlap = 40)
save<-saveRDS(D80milo_k45, file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D80_k45_pca10_d10_prop0.1_2026_04_17.rds")

#Run DA glmm
##run tasks in a serial manner
bpparam <- SerialParam()
register(bpparam)
##use non-negative least squares (NNLS) Haseman-Elston solver,REML is Restricted Maximum Likelihood),fdr.weighting=spatial FDR weighting scheme to use
print("D80 k75")
glm_results_D80 <- testNhoods(D80milo_k75, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)
which(checkSeparation(D80milo_k75, design.df=pbmc.design, condition="time_point", min.val=10))


da_results_D80 <- testNhoods(D80milo_k75, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)



comp.da <- merge(da_results_D80, glm_results_D80, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()



## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_NhoodGraph_2026_04_17.pdf",width=10, height=10)
plotNhoodGraphDA(D80milo_k75, da_results_D80, layout="UMAP",alpha=0.1) 
dev.off()

# Add cell types to results 
da_results_D80 <- annotateNhoods(D80milo_k75, da_results_D80, coldata_col = "celltypes_greek")
head(da_results_D80)

saveRDS(da_results_D80, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_2026_04_17.rds")
da_results_D80<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_2026_04_17.rds")
write.table(da_results_D80, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_2026_04_17.txt", sep="\t", row.names=FALSE)
rm(comp.da)
### Try with K45
print("D80 k45")

which(checkSeparation(D80milo_k45, design.df=pbmc.design, condition="time_point", min.val=10))

glm_results_D80_k45 <- testNhoods(D80milo_k45, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)

da_results_D80_k45 <- testNhoods(D80milo_k45, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)

comp.da <- merge(da_results_D80_k45, glm_results_D80_k45, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D80_k45_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()

## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D80_k45_DA_glmm_NhoodGraph_2026_04_17.pdf")
plotNhoodGraphDA(D80milo_k45,da_results_D80_k45, layout="UMAP",alpha=0.1) 
dev.off()

# Add cell types to results 
da_results_D80_k45 <- annotateNhoods(D80milo_k45, da_results_D80_k45, coldata_col = "celltypes_greek")
head(da_results_D80_k45)

saveRDS(da_results_D80_k45, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k45_DA_glmm_2026_04_17.rds")


print("Run D8 vs D2")
#Split by time points
sce_D82 <- sce[, sce$time_point %in% c("D2", "D8")]

#Make milo object from the SCE object
D82milo <- Milo(sce_D82) 


#Add UMAP from SCE object to milo object
reducedDim(D82milo, "UMAP") <- reducedDim(sce_D82, "UMAP")


####=======  k=# of nearest neighbors to include(# of samples x 5 = 75) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 5=75), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D82milo_k75 <- buildGraph(D82milo, k = 75, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D82milo_k75 <- makeNhoods(D82milo_k75, prop = 0.1, k = 75, d = 10, refined = TRUE, reduced_dims = "UMAP",refinement_scheme="graph") 

#add animal_time to metadata
colData(D82milo_k75)$ObsID <- paste(colData(D82milo_k75)$animal, colData(D82milo_k75)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D82milo_k75 <- countCells(D82milo_k75, meta.data = data.frame(colData(D82milo_k75)), samples="Animal_time") 
D82milo_k75

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k75Ncounts<-nhoodCounts(D82milo_k75)
k75Ncounts_s<-saveRDS(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D8v2_2026_04_17.rds")
k75Ncountsm<-as.matrix(k75Ncounts)
k75Ncountsdf<-as.data.frame(k75Ncountsm)
k75Ncounts<-k75Ncountsdf
k75Ncounts_s<-write.csv(k75Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k75Ncounts_D8v2_2026_04_17.csv")
avg_k75Ncounts<-colMeans(k75Ncounts)
print(avg_k75Ncounts)
median_k75Ncounts<-sapply(k75Ncounts,median)
print(median_k75Ncounts)
IQR_k75Ncounts<-sapply(k75Ncounts,IQR)
print(IQR_k75Ncounts)
save<-saveRDS(D82milo_k75,file="/scRNAseq/Sal_5pigs_2026/DA/milo_D8v2_k75_pca10_d10_prop0.1_2026_04_17.rds")
####=======  k=# of nearest neighbors to include(# of samples x 3 = 45) =======####
#Make K-nearest Neighbor Graph, k=# of nearest neighbors to include(s# x 3=45), d= # of dimensions/pca to use (use the same pca as was used in Seurat)
D82milo_k45 <- buildGraph(D82milo, k = 45, d = 10)

#Define neighborhoods on UMAP, prop=proportion of graph vertices to randomly sample, k & d need to be the same as buildGraph
D82milo_k45 <- makeNhoods(D82milo_k45, prop = 0.1, k = 45, d = 10, refined = TRUE, reduced_dims = "UMAP",refinement_scheme="graph") 

#add animal_time to metadata
colData(D82milo_k45)$ObsID <- paste(colData(D82milo_k45)$animal, colData(D82milo_k45)$time_point, sep="_")

#Make the nhood X sample counts matrix, stored in milo object 
D82milo_k45 <- countCells(D82milo_k45, meta.data = data.frame(colData(D82milo_k45)), samples="Animal_time") 
D82milo_k45

#generate table to ensure there are 3-5 cells from each sample in each neighborhood
k45Ncounts<-nhoodCounts(D82milo_k45)
k45Ncounts_s<-saveRDS(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/k45Ncounts_D8v2_2026_04_17.rds")
k45Ncountsm<-as.matrix(k45Ncounts)
k45Ncountsdf<-as.data.frame(k45Ncountsm)
k45Ncounts<-k45Ncountsdf
k45Ncounts_s<-write.csv(k45Ncounts,file="/scRNAseq/Sal_5pigs_2026/DA/k45/k45Ncounts_D8v2_2026_04_17.csv")
avg_k45Ncounts<-colMeans(k45Ncounts)
print(avg_k45Ncounts)
median_k45Ncounts<-sapply(k45Ncounts,median)
print(median_k45Ncounts)
IQR_k45Ncounts<-sapply(k45Ncounts,IQR)
print(IQR_k45Ncounts)

det_K<-rbind(avg_k75Ncounts,median_k75Ncounts,IQR_k75Ncounts,avg_k45Ncounts,median_k45Ncounts,IQR_k45Ncounts)

det_k_s<-write.csv(det_K, file="/scRNAseq/Sal_5pigs_2026/DA/K45_75_D8v2_avg_IQR_2026_04_17.csv")

k45Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D82milo_k45)),
                        "Var"=rowVars(nhoodCounts(D82milo_k45)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D82_k45_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k45Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

k75Ncounts.var.df <- data.frame("Mean"=rowMeans2(nhoodCounts(D82milo_k75)),
                        "Var"=rowVars(nhoodCounts(D82milo_k75)))
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_Var_vs_Mean_2026_04_17.pdf")
ggplot(k75Ncounts.var.df, aes(x=Mean, y=Var)) +
    geom_point() +
    geom_smooth() +
    geom_abline(lty=2, colour='red') +
    theme_cowplot() +
    NULL
dev.off()

#note by Carrie: counts= # of neighbourhoods neighbourhood size = # of cells in neighbourhood
pdf("/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D82milo_k75)
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D82_k45_pca10_d10_prop0.1_plotNhoodSizeHist_2026_04_17.pdf")
miloR::plotNhoodSizeHist(D82milo_k45)
dev.off()


#Design model matrix from metadata, make sure the rownames match with nhoodCounts column names & are in same order
pbmc.design <- data.frame(colData(D82milo_k75))[,c("Animal_time", "time_point", "Animal")]
pbmc.design$Animal <- as.factor(pbmc.design$Animal) 
pbmc.design <- distinct(pbmc.design)
rownames(pbmc.design) <- pbmc.design$Animal_time
pbmc.design <- pbmc.design[colnames(nhoodCounts(D82milo_k75)), , drop=FALSE]
head(pbmc.design)

#Calculate Nhood distance to later use to account for the overlap between neighbourhoods & build Nhood graph
D82milo_k75 <- calcNhoodDistance(D82milo_k75, d = 10)
D82milo_k75 <- buildNhoodGraph(D82milo_k75, overlap=70)
save<-saveRDS(D82milo_k75, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_2026_04_17.rds")
D82milo_k75 <- buildNhoodGraph(D82milo_k75, overlap=60)
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_NhoodGraph_2026_04_17_overlap60.pdf",width=10, height=10)
plotNhoodGraph(D82milo_k75, layout="UMAP")
dev.off()
D82milo_k75 <- buildNhoodGraph(D82milo_k75, overlap=50)
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_NhoodGraph_2026_04_17_overlap50.pdf",width=10, height=10)
plotNhoodGraph(D82milo_k75, layout="UMAP")
dev.off()
D82milo_k75 <- buildNhoodGraph(D82milo_k75, overlap=70)
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_NhoodGraph_2026_04_17_overlap70.pdf",width=10, height=10)
plotNhoodGraph(D82milo_k75, layout="UMAP")
dev.off()


D82milo_k45 <- calcNhoodDistance(D82milo_k45, d = 10)
D82milo_k45 <- buildNhoodGraph(D82milo_k45,overlap=40)
save<-saveRDS(D82milo_k45, file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D82_k45_pca10_d10_prop0.1_2026_04_17.rds")


#Run DA glmm
##run tasks in a serial manner
bpparam <- SerialParam()
register(bpparam)
##use non-negative least squares (NNLS) Haseman-Elston solver,REML is Restricted Maximum Likelihood),fdr.weighting=spatial FDR weighting scheme to use
print("D82 k75")
glm_results_D82 <- testNhoods(D82milo_k75, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)
which(checkSeparation(D82milo_k75, design.df=pbmc.design, condition="time_point", min.val=10))


da_results_D82 <- testNhoods(D82milo_k75, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)



comp.da <- merge(da_results_D82, glm_results_D82, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()



## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_NhoodGraph_2026_04_17.pdf",width=10, height=10)
plotNhoodGraphDA(D82milo_k75, da_results_D82, layout="UMAP",alpha=0.1) 
dev.off()

# Add cell types to results 
da_results_D82 <- annotateNhoods(D82milo_k75, da_results_D82, coldata_col = "celltypes_greek")
head(da_results_D82)

saveRDS(da_results_D82, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_2026_04_17.rds")
write.table(da_results_D82, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_2026_04_17.txt", sep="\t", row.names=FALSE)


rm(comp.da)
### Try with K45
print("D82 k45")

which(checkSeparation(D82milo_k45, design.df=pbmc.design, condition="time_point", min.val=10))

glm_results_D82_k45 <- testNhoods(D82milo_k45, design = ~ time_point + Animal, design.df = pbmc.design, fdr.weighting="graph-overlap",REML=TRUE, norm.method="TMM", BPPARAM = bpparam)

da_results_D82_k45 <- testNhoods(D82milo_k45, design = ~ time_point + (1|Animal), design.df = pbmc.design, fdr.weighting="graph-overlap",glmm.solver="HE-NNLS", REML=TRUE, norm.method="TMM", BPPARAM = bpparam, fail.on.error=FALSE, force=TRUE)

comp.da <- merge(da_results_D82_k45, glm_results_D82_k45, by='Nhood')
comp.da$Sig <- "none"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y < 0.1] <- "Both"
comp.da$Sig[comp.da$SpatialFDR.x >= 0.1 & comp.da$SpatialFDR.y < 0.1] <- "GLM"
comp.da$Sig[comp.da$SpatialFDR.x < 0.1 & comp.da$SpatialFDR.y >= 0.1] <- "GLMM"

pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D82_k45_glm_vs_glmm_2026_04_17.pdf")
ggplot(comp.da, aes(x=logFC.x, y=logFC.y)) +geom_point(data=comp.da[, c("logFC.x", "logFC.y")], aes(x=logFC.x, y=logFC.y),colour='grey80', size=1) + geom_point(aes(colour=Sig)) +labs(x="GLMM LFC", y="GLM LFC") +facet_wrap(~Sig) 
dev.off()

## Plot DA neighbourhood graph
pdf(file="/scRNAseq/Sal_5pigs_2026/DA/k45/milo_D82_45_DA_glmm_NhoodGraph_2026_04_17.pdf")
plotNhoodGraphDA(D82milo_k45,da_results_D82_k45, layout="UMAP",alpha=0.1) 
dev.off()

# Add cell types to results 
da_results_D82_k45 <- annotateNhoods(D82milo_k45, da_results_D82_k45, coldata_col = "celltypes_greek")
head(da_results_D82_k45)

saveRDS(da_results_D82_k45, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k45_DA_glmm_2026_04_17.rds")

sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux 9.6 (Plow)

# Matrix products: default
# BLAS/LAPACK: /micromamba/envs/miloR2.0/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0

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
#  [1] Cairo_1.6-2                 cowplot_1.1.3              
#  [3] BiocParallel_1.36.0         irlba_2.3.5.1              
#  [5] Matrix_1.6-5                scRNAseq_2.16.0            
#  [7] patchwork_1.3.0             dplyr_1.1.4                
#  [9] scran_1.30.2                scater_1.30.1              
# [11] ggplot2_3.5.1               scuttle_1.12.0             
# [13] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [15] Biobase_2.62.0              GenomicRanges_1.54.1       
# [17] GenomeInfoDb_1.38.8         IRanges_2.36.0             
# [19] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# [21] MatrixGenerics_1.14.0       matrixStats_1.4.1          
# [23] miloR_2.1.3                 edgeR_4.0.16               
# [25] limma_3.58.1               

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3            magrittr_2.0.3               
#   [3] GenomicFeatures_1.54.4        ggbeeswarm_0.7.2             
#   [5] farver_2.1.2                  BiocIO_1.12.0                
#   [7] zlibbioc_1.48.2               vctrs_0.6.5                  
#   [9] Rsamtools_2.18.0              memoise_2.0.1                
#  [11] DelayedMatrixStats_1.24.0     RCurl_1.98-1.16              
#  [13] progress_1.2.3                htmltools_0.5.8.1            
#  [15] S4Arrays_1.2.1                AnnotationHub_3.10.1         
#  [17] curl_6.2.2                    BiocNeighbors_1.20.2         
#  [19] SparseArray_1.2.4             pracma_2.4.4                 
#  [21] cachem_1.1.0                  GenomicAlignments_1.38.2     
#  [23] igraph_2.1.2                  mime_0.12                    
#  [25] lifecycle_1.0.4               pkgconfig_2.0.3              
#  [27] rsvd_1.0.5                    R6_2.5.1                     
#  [29] fastmap_1.2.0                 GenomeInfoDbData_1.2.11      
#  [31] shiny_1.8.1.1                 digest_0.6.36                
#  [33] numDeriv_2016.8-1.1           colorspace_2.1-1             
#  [35] AnnotationDbi_1.64.1          dqrng_0.4.1                  
#  [37] ExperimentHub_2.10.0          RSQLite_2.3.7                
#  [39] beachmat_2.18.1               labeling_0.4.3               
#  [41] filelock_1.0.3                mgcv_1.9-1                   
#  [43] httr_1.4.7                    polyclip_1.10-7              
#  [45] abind_1.4-8                   compiler_4.3.3               
#  [47] bit64_4.0.5                   withr_3.0.2                  
#  [49] viridis_0.6.5                 DBI_1.2.3                    
#  [51] ggforce_0.4.2                 biomaRt_2.58.2               
#  [53] MASS_7.3-60                   rappdirs_0.3.3               
#  [55] DelayedArray_0.28.0           rjson_0.2.21                 
#  [57] bluster_1.12.0                gtools_3.9.5                 
#  [59] tools_4.3.3                   vipor_0.4.7                  
#  [61] beeswarm_0.4.0                interactiveDisplayBase_1.40.0
#  [63] httpuv_1.6.15                 glue_1.8.0                   
#  [65] restfulr_0.0.15               nlme_3.1-165                 
#  [67] promises_1.3.0                grid_4.3.3                   
#  [69] cluster_2.1.6                 generics_0.1.3               
#  [71] gtable_0.3.6                  ensembldb_2.26.0             
#  [73] tidyr_1.3.1                   hms_1.1.3                    
#  [75] xml2_1.3.6                    BiocSingular_1.18.0          
#  [77] tidygraph_1.3.1               ScaledMatrix_1.10.0          
#  [79] metapod_1.10.1                XVector_0.42.0               
#  [81] ggrepel_0.9.6                 BiocVersion_3.18.1           
#  [83] pillar_1.10.0                 stringr_1.5.1                
#  [85] later_1.3.2                   splines_4.3.3                
#  [87] tweenr_2.0.3                  BiocFileCache_2.10.2         
#  [89] lattice_0.22-6                rtracklayer_1.62.0           
#  [91] bit_4.0.5                     tidyselect_1.2.1             
#  [93] locfit_1.5-9.10               Biostrings_2.70.3            
#  [95] gridExtra_2.3                 ProtGenerics_1.34.0          
#  [97] graphlayouts_1.2.1            statmod_1.5.0                
#  [99] stringi_1.8.4                 lazyeval_0.2.2               
# [101] yaml_2.3.9                    codetools_0.2-20             
# [103] ggraph_2.2.1                  tibble_3.2.1                 
# [105] BiocManager_1.30.23           cli_3.6.3                    
# [107] xtable_1.8-4                  munsell_0.5.1                
# [109] Rcpp_1.0.13-1                 dbplyr_2.5.0                 
# [111] png_0.1-8                     XML_3.99-0.17                
# [113] parallel_4.3.3                blob_1.2.4                   
# [115] prettyunits_1.2.0             AnnotationFilter_1.26.0      
# [117] sparseMatrixStats_1.14.0      bitops_1.0-9                 
# [119] viridisLite_0.4.2             scales_1.3.0                 
# [121] purrr_1.0.2                   crayon_1.5.3                 
# [123] rlang_1.1.4                   KEGGREST_1.42.0              