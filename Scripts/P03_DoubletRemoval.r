.libPaths()
#[1] "/micromamba/envs/SCE/lib/R/library"
{set.seed(123)
library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)
library(ggplot2)
library(BiocParallel)}



All_sce <- readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_SCE_QC_2026_04_15.rds")
All_sce 

#Find Doublets in each individual sample
All_sce.2 <- scDblFinder(All_sce, multiSampleMode ="split",samples="SampleID", BPPARAM=MulticoreParam(3))
#Check distrublution of doublets??
table(All_sce.2$scDblFinder.class)
#singlet doublet 
# 151808   30407 
print("table(All_sce.2$scDblFinder.sample)")
table(All_sce.2$scDblFinder.sample)
#S01_852_D0_S S02_852_D2_S S03_852_D8_S S04_854_D0_S S05_854_D2_S S06_854_D8_S 
#   10964         5069         6480         1322         7320         3729 
#S07_864_D0_S S08_864_D2_S S09_864_D8_S S10_842_D0_S S11_842_D2_S S12_842_D8_S 
#    3747        13167        10993        25649        14816        24856 
#S13_853_D0_S S14_853_D2_S S15_853_D8_S 
#   33842         6116        14145 #Save doublet analysis
All_sce.2
# dim: 22045 61884
saveRDS(All_sce.2,file ="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_QC_IDed_Doublets_sce_2026_04_15.rds")


class_counts <- table(All_sce.2$scDblFinder.class)
pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/5pigs_Distribution_Singlets_Doublets.pdf")
barplot(class_counts, main = "Distribution of Singlets and Doublets", xlab = "Class", ylab = "Count", col = c("green", "red"))
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/5pigs_Histogram_Doublet_Scores.pdf")
ggplot(data.frame(Score = All_sce.2$scDblFinder.score), aes(x = Score)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  labs(title = "Histogram of Doublet Scores", x = "Doublet Score", y = "Frequency")
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/5pigs_Histogram_Weighted_Doublet_Score.pdf")
ggplot(data.frame(WeightedScore = All_sce.2$scDblFinder.weighted), aes(x = WeightedScore)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  labs(title = "Histogram of Weighted Scores", x = "Weighted Score", y = "Frequency")
dev.off()


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
#  [1] BiocParallel_1.36.0         scDblFinder_1.16.0
#  [3] tidyr_1.3.1                 PCAtools_2.14.0
#  [5] ggrepel_0.9.6               scater_1.30.1
#  [7] ggplot2_3.5.1               Matrix_1.6-5
#  [9] scran_1.30.0                scuttle_1.12.0
# [11] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [13] Biobase_2.62.0              GenomicRanges_1.54.1
# [15] GenomeInfoDb_1.38.1         IRanges_2.36.0
# [17] S4Vectors_0.40.2            BiocGenerics_0.48.1
# [19] MatrixGenerics_1.14.0       matrixStats_1.5.0

# loaded via a namespace (and not attached):
#  [1] bitops_1.0-9              gridExtra_2.3
#  [3] rlang_1.1.5               magrittr_2.0.3
#  [5] compiler_4.3.3            DelayedMatrixStats_1.24.0
#  [7] vctrs_0.6.5               reshape2_1.4.4
#  [9] stringr_1.5.1             pkgconfig_2.0.3
# [11] crayon_1.5.3              XVector_0.42.0
# [13] labeling_0.4.3            Rsamtools_2.18.0
# [15] ggbeeswarm_0.7.2          purrr_1.0.2
# [17] bluster_1.12.0            zlibbioc_1.48.0
# [19] beachmat_2.18.0           jsonlite_1.8.9
# [21] DelayedArray_0.28.0       irlba_2.3.5.1
# [23] parallel_4.3.3            cluster_2.1.8
# [25] R6_2.5.1                  stringi_1.8.4
# [27] limma_3.58.1              rtracklayer_1.62.0
# [29] xgboost_2.1.4.1           Rcpp_1.0.14
# [31] igraph_2.0.3              tidyselect_1.2.1
# [33] abind_1.4-5               yaml_2.3.10
# [35] viridis_0.6.5             codetools_0.2-20
# [37] lattice_0.22-6            tibble_3.2.1
# [39] plyr_1.8.9                withr_3.0.2
# [41] Biostrings_2.70.1         pillar_1.10.1
# [43] generics_0.1.3            RCurl_1.98-1.16
# [45] sparseMatrixStats_1.14.0  munsell_0.5.1
# [47] scales_1.3.0              glue_1.8.0
# [49] metapod_1.10.0            tools_4.3.3
# [51] BiocIO_1.12.0             BiocNeighbors_1.20.0
# [53] data.table_1.16.4         ScaledMatrix_1.10.0
# [55] locfit_1.5-9.11           GenomicAlignments_1.38.0
# [57] XML_3.99-0.17             cowplot_1.1.3
# [59] grid_4.3.3                edgeR_4.0.16
# [61] colorspace_2.1-1          GenomeInfoDbData_1.2.11
# [63] beeswarm_0.4.0            BiocSingular_1.18.0
# [65] restfulr_0.0.15           vipor_0.4.7
# [67] cli_3.6.3                 rsvd_1.0.5
# [69] S4Arrays_1.2.0            viridisLite_0.4.2
# [71] dplyr_1.1.4               gtable_0.3.6
# [73] SparseArray_1.2.2         dqrng_0.3.2
# [75] rjson_0.2.23              farver_2.1.2
# [77] lifecycle_1.0.4           statmod_1.5.0
# [79] MASS_7.3-60.0.1