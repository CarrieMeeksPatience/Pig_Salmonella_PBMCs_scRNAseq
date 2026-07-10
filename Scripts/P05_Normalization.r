#Sources for code:
## https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html
##https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/06_FeatureSelectionAndDimensionalityReduction.html#dimensionality-reduction
.libPaths()
#[1] "/micromamba/envs/SCE/lib/R/library"

{set.seed(123)
library(scuttle)
library(scran)
library(SingleCellExperiment)
library(Matrix)
library(scater)
library(PCAtools)
library(tidyr)}

sce<-readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_SCE_2026_04_15.rds")


#filter out any remaining poor quality cells via scuttle
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
sce <- sce[,!qcfilter$discard]
summary(qcfilter$discard)
 #Mode     FALSE 
 #logical  84569
#Cluster similar cells based on their expression profiles using either log-expression values via scran
clusters <- quickCluster(sce)
#compute scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools via scran
sce <- computeSumFactors(sce, clusters=clusters)
#Gets the size factors for all cells via SingleCellExperiment
summary(sizeFactors(sce))
### Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
## 0.2785  0.5464  0.7753  1.0000  1.2035 14.4291
#Compute log-transformed normalized expression values via scuttle
sce <- logNormCounts(sce)

#Check variance of the data
#Model the variance of the log-expression profiles for each gene, decomposing it into technical and biological components based on a fitted mean-variance trend via scran
dec <- modelGeneVar(sce)

pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/Scran_Model_Gene_Var.pdf")
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
dev.off()

# Get the HVGs for clustering
#Get the top 10% of genes.
top.hvgs <- getTopHVGs(dec, prop=0.1)
print("length(top.hvgs),prop=0.1")
length(top.hvgs)
#[1] 1280 genes(lower than 10% of genes becuase Genes with very low expression across cells may be excluded from the analysis, as they are unlikely to be informative for downstream analyses.)
top.hvgs_s<-dput(top.hvgs, file="/scRNAseq/Sal_5pigs_2026/Normalization/HVGs_scran.R")
##Plot the HVGs
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/Scran_top1280HVGs_vplot.pdf")
plotExpression(sce, features = top.hvgs[1:20], point_alpha = 0.05, jitter = "jitter")
dev.off()
# Perform PCA on the HVGs
sce <- fixedPCA(sce, subset.row=top.hvgs) #50 PCs

## ensure PCA is there
reducedDimNames(sce)
## [1] "PCA"
## check the dimensions of the PCA matrix
dim(reducedDim(sce, "PCA"))
#[1] 84569     50

## check the variance explained by each PC
percent.var <- attr(reducedDim(sce), "percentVar")
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_elbow_plot.pdf")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
dev.off()
sce_s<-saveRDS(sce,file="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_ScranNorm_SCE_2026_04_15.rds")

## make elbow plot with 0.5%, 0.4%, 0.3%, 0.25% & 0.2%  of the variance marked
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/elbow_plot_varianceMarked.pdf")
# Generate the plot with custom x-axis
plot(percent.var, log = "y", xlab = "PC", ylab = "Variance explained (%)", xaxt = 'n')

# Add horizontal lines at 0.5%, 0.4%, 0.3%, 0.25%, and 0.2% variance
abline(h = 0.5, col = "blue", lty = 2)
abline(h = 0.45, col = "orange", lty = 2) 
abline(h = 0.4, col = "cyan", lty = 2)
abline(h = 0.3, col = "purple", lty = 2)
abline(h = 0.25, col = "green", lty = 2)
abline(h = 0.2, col = "red", lty = 2)

# Customize the x-axis ticks and labels to be in increments of 5
axis(1, at = seq(0, length(percent.var), by = 5))

# Add a legend to the plot
legend("topright", legend = c("0.5%", "0.45%","0.4%", "0.3%", "0.25%", "0.2%"), 
       col = c("blue", "orange","cyan", "purple", "green", "red"), lty = 2)
dev.off()

# Make a data frame of the variance explained by each PC
    #Calculate the cumulative percentage of variance explained
cumu <- cumsum(percent.var)
co1 <- which(cumu > 90 & percent.var < 5)[1]
co2 <- sort(which((percent.var[1:length(percent.var) - 1] - percent.var[2:length(percent.var)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(percent.var = percent.var, cumu = cumu, rank = 1:length(percent.var))
plot_df 
write.csv(plot_df, file = "/scRNAseq/Sal_5pigs_2026/Normalization/scran_variance_explained.csv")
# Plot the variance explained by each PC
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_variance_explained.pdf",width=15, height=10)
ggplot(plot_df, aes(cumu, percent.var, label = rank, color = rank > pcs)) +
  geom_text() + geom_vline(xintercept = 33, color = "blue") +
  geom_hline(yintercept = min(percent.var[percent.var > 0.4]), color = "red") + theme_bw()+scale_color_manual(values = c("33% Cumulative Variance" = "blue", "Min Variance > 0.4%" = "red", "TRUE" = "black", "FALSE" = "black")) +
  labs(color = "Legend") +
  theme(legend.position = "topright") + xlab("Cumulative Variance (%)") + ylab("Variance Explained (%)")
dev.off()

## plot PCs
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_PCs_plot.pdf")
plotReducedDim(sce, dimred="PCA", ncomponents=6,
    colour_by="SampleID")
dev.off()

## Plot correlations between different variables and our PC scores
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_PCs_correlations.pdf",width=15, height=10)
explan_pcs <- getExplanatoryPCs(sce,n_dimred=50,
    variables = c(
        "nCount_RNA",
        "nFeature_RNA",
        "SampleID",
        "Animal",
        "percent_mito_ribo"
    )
)
plotExplanatoryPCs(explan_pcs/100)
dev.off()
# distribution of correlations between each gene's expression and our variables of interest
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_PCs_density.pdf", width=20, height=10)
plotExplanatoryVariables(sce,
                         variables = c(
                         "nCount_RNA",
                        "nFeature_RNA",
                        "SampleID",
                        "Animal",
                        "percent_mito_ribo"
                                        ))
plotExplanatoryPCs(explan_pcs/100)
dev.off()


# extract variance explained for PCAtools
pca_pct_variance <- data.frame(variance = attr(reducedDim(sce, "PCA"), "percentVar"))
pca_pct_variance$PC <- 1:nrow(pca_pct_variance)
#How many PCs explain more than 1% of the variance in the entire dataset.
table(pca_pct_variance$variance > 1)
##FALSE  TRUE
##  45     5
#So PC 4 explain more than 1% of the variance in the entire dataset
# Choose # of PCs to use via PCAtools
chosen_elbow <- findElbowPoint(pca_pct_variance$variance)
chosen_elbow
# [1] 4 ; 4 PCs are reccomended to be informative
#plot PCAtools elbow
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scran_PCAtools_elbow.pdf")
pca_pct_variance %>% ggplot(aes(PC, variance)) + geom_point() +geom_vline(xintercept = chosen_elbow)
dev.off()
## run denoise PCA step to Denoise log-expression data by removing principal components corresponding to technical noise. This gives a rough idea of the number of PCs that are likely to be informative.
sce2 <- denoisePCA(sce, technical = dec)
# check dimensions of the "denoised" PCA
ncol(reducedDim(sce2, "PCA"))
# [1] 5 ; 5 PCs are reccomended to be informative
#OR denoise till 1% of varince is explained
sce3 <- denoisePCA(sce, dec, subset.row=top.hvgs)
ncol(reducedDim(sce3, "PCA"))
# [1] 5 ; 5 PCs are reccomended to be informative

sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux 9.6 (Plow)


# Matrix products: default
# BLAS/LAPACK: /micromamba/envs/SCE/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0

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
#  [1] tidyr_1.3.1                 PCAtools_2.14.0
#  [3] ggrepel_0.9.6               scater_1.30.1
#  [5] ggplot2_3.5.1               Matrix_1.6-5
#  [7] scran_1.30.0                scuttle_1.12.0
#  [9] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [11] Biobase_2.62.0              GenomicRanges_1.54.1
# [13] GenomeInfoDb_1.38.1         IRanges_2.36.0
# [15] S4Vectors_0.40.2            BiocGenerics_0.48.1
# [17] MatrixGenerics_1.14.0       matrixStats_1.5.0

# loaded via a namespace (and not attached):
#  [1] beeswarm_0.4.0            gtable_0.3.6
#  [3] lattice_0.22-6            vctrs_0.6.5
#  [5] tools_4.3.3               bitops_1.0-9
#  [7] generics_0.1.3            parallel_4.3.3
#  [9] tibble_3.2.1              cluster_2.1.8
# [11] pkgconfig_2.0.3           BiocNeighbors_1.20.0
# [13] sparseMatrixStats_1.14.0  dqrng_0.3.2
# [15] lifecycle_1.0.4           GenomeInfoDbData_1.2.11
# [17] farver_2.1.2              stringr_1.5.1
# [19] compiler_4.3.3            statmod_1.5.0
# [21] munsell_0.5.1             bluster_1.12.0
# [23] codetools_0.2-20          vipor_0.4.7
# [25] RCurl_1.98-1.16           pillar_1.10.1
# [27] crayon_1.5.3              BiocParallel_1.36.0
# [29] DelayedArray_0.28.0       limma_3.58.1
# [31] viridis_0.6.5             abind_1.4-5
# [33] tidyselect_1.2.1          rsvd_1.0.5
# [35] locfit_1.5-9.11           metapod_1.10.0
# [37] stringi_1.8.4             purrr_1.0.2
# [39] reshape2_1.4.4            BiocSingular_1.18.0
# [41] dplyr_1.1.4               labeling_0.4.3
# [43] cowplot_1.1.3             grid_4.3.3
# [45] colorspace_2.1-1          cli_3.6.3
# [47] SparseArray_1.2.2         magrittr_2.0.3
# [49] S4Arrays_1.2.0            edgeR_4.0.16
# [51] withr_3.0.2               DelayedMatrixStats_1.24.0
# [53] scales_1.3.0              ggbeeswarm_0.7.2
# [55] XVector_0.42.0            igraph_2.0.3
# [57] gridExtra_2.3             ScaledMatrix_1.10.0
# [59] beachmat_2.18.0           viridisLite_0.4.2
# [61] irlba_2.3.5.1             rlang_1.1.5
# [63] Rcpp_1.0.14               glue_1.8.0
# [65] plyr_1.8.9                R6_2.5.1
# [67] zlibbioc_1.48.0
