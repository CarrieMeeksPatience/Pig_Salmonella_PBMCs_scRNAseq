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


library(ggrepel)

pca_data <- Embeddings(All3_integrated, "pca")[, 1:2] |>
  as.data.frame()

pca_data$Animal <- All3_integrated$Animal
pca_data$DPI <- All3_integrated$DPI

# One point per Animal x DPI
pca_summary <- pca_data |>
  group_by(Animal, DPI) |>
  summarise(
    PC_1 = mean(PC_1, na.rm = TRUE),
    PC_2 = mean(PC_2, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

pdf("/scRNAseq/Sal_5pigs_2026/DE/PCA_AnimalDPI_seurat.pdf", width = 6, height = 6)
    ggplot(pca_summary, aes(
    x = PC_1,
    y = PC_2,
    color = Animal,
    shape = DPI
  )) +
    geom_point(size = 3) +
stat_ellipse(aes(group = DPI, color = DPI),type = "t",linetype = 2,level = 0.68)+
 coord_fixed()+
    theme_classic()+
    theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Align title to the middle
    plot.title.position = "plot")+
    labs(
      title = "PCA of scRNAseq Data",
      x = "PC_1",
      y = "PC_2",
      shape = "DPI",
      colour = "Pig ID"
    )
dev.off()

sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: AlmaLinux 9.6 (Sage Margay)

# Matrix products: default
# BLAS/LAPACK: /micromamba/envs/seurat+milo.new/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0

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
#  [1] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
#  [3] Biobase_2.62.0              GenomicRanges_1.54.1       
#  [5] GenomeInfoDb_1.38.8         IRanges_2.36.0             
#  [7] S4Vectors_0.40.2            BiocGenerics_0.48.1        
#  [9] MatrixGenerics_1.14.0       matrixStats_1.5.0          
# [11] findPC_1.0                  PCAtools_2.14.0            
# [13] ggrepel_0.9.6               RColorBrewer_1.1-3         
# [15] clustree_0.5.1              ggraph_2.2.2               
# [17] scCustomize_3.3.0           tidyr_1.3.2                
# [19] scIntegrationMetrics_1.2.0  STACAS_2.4.1               
# [21] dplyr_1.1.4                 ggplot2_4.0.3              
# [23] Seurat_5.5.0                SeuratObject_5.4.0         
# [25] sp_2.2-0                   

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23          splines_4.3.3            
#   [3] later_1.4.4               bitops_1.1-0             
#   [5] tibble_3.3.0              R.oo_1.27.1              
#   [7] polyclip_1.10-7           janitor_2.2.1            
#   [9] fastDummies_1.7.6         lifecycle_1.0.5          
#  [11] globals_0.19.1            lattice_0.22-7           
#  [13] MASS_7.3-60               magrittr_2.0.4           
#  [15] plotly_4.12.1             httpuv_1.6.16            
#  [17] sctransform_0.4.2         spam_2.11-1              
#  [19] spatstat.sparse_3.1-0     reticulate_1.43.0        
#  [21] cowplot_1.2.0             pbapply_1.7-4            
#  [23] lubridate_1.9.5           zlibbioc_1.48.2          
#  [25] abind_1.4-8               Rtsne_0.17               
#  [27] purrr_1.1.0               R.utils_2.13.0           
#  [29] RCurl_1.98-1.19           tweenr_2.0.3             
#  [31] GenomeInfoDbData_1.2.11   circlize_0.4.18          
#  [33] irlba_2.3.7               listenv_1.0.0            
#  [35] spatstat.utils_3.2-0      vegan_2.7-5              
#  [37] goftest_1.2-3             RSpectra_0.16-2          
#  [39] dqrng_0.4.1               spatstat.random_3.4-2    
#  [41] fitdistrplus_1.2-6        parallelly_1.45.1        
#  [43] DelayedMatrixStats_1.24.0 permute_0.9-10           
#  [45] codetools_0.2-20          DelayedArray_0.28.0      
#  [47] ggforce_0.5.0             tidyselect_1.2.1         
#  [49] shape_1.4.6.1             farver_2.1.2             
#  [51] ScaledMatrix_1.10.0       viridis_0.6.5            
#  [53] spatstat.explore_3.5-3    jsonlite_2.0.0           
#  [55] BiocNeighbors_1.20.2      tidygraph_1.3.1          
#  [57] progressr_1.0.0           ggridges_0.5.7           
#  [59] survival_3.8-3            tools_4.3.3              
#  [61] ica_1.0-3                 Rcpp_1.1.0               
#  [63] glue_1.8.0                gridExtra_2.3.1          
#  [65] SparseArray_1.2.4         mgcv_1.9-3               
#  [67] withr_3.0.3               fastmap_1.2.0            
#  [69] mcprogress_0.1.1          rsvd_1.0.5               
#  [71] digest_0.6.37             timechange_0.4.0         
#  [73] R6_2.6.1                  mime_0.13                
#  [75] ggprism_1.0.7             colorspace_2.1-2         
#  [77] scattermore_1.2           tensor_1.5.1             
#  [79] spatstat.data_3.1-9       R.methodsS3_1.8.2        
#  [81] generics_0.1.4            data.table_1.17.8        
#  [83] graphlayouts_1.2.5        httr_1.4.8               
#  [85] htmlwidgets_1.6.4         S4Arrays_1.2.1           
#  [87] uwot_0.2.4                pkgconfig_2.0.3          
#  [89] gtable_0.3.6              lmtest_0.9-40            
#  [91] S7_0.2.2                  XVector_0.42.0           
#  [93] htmltools_0.5.9           dotCall64_1.2            
#  [95] scales_1.4.0              png_0.1-9                
#  [97] spatstat.univar_3.1-4     snakecase_0.11.1         
#  [99] reshape2_1.4.5            nlme_3.1-168             
# [101] zoo_1.8-14                cachem_1.1.0             
# [103] GlobalOptions_0.1.4       stringr_1.6.0            
# [105] KernSmooth_2.23-26        parallel_4.3.3           
# [107] miniUI_0.1.2              vipor_0.4.7              
# [109] ggrastr_1.0.2             pillar_1.11.1            
# [111] grid_4.3.3                vctrs_0.6.5              
# [113] RANN_2.6.2                promises_1.3.3           
# [115] BiocSingular_1.18.0       beachmat_2.18.1          
# [117] xtable_1.8-8              cluster_2.1.8.1          
# [119] beeswarm_0.4.0            paletteer_1.7.0          
# [121] cli_3.6.5                 compiler_4.3.3           
# [123] rlang_1.1.6               crayon_1.5.3             
# [125] future.apply_1.20.2       labeling_0.4.3           
# [127] rematch2_2.1.2            plyr_1.8.9               
# [129] forcats_1.0.1             ggbeeswarm_0.7.3         
# [131] stringi_1.8.7             viridisLite_0.4.3        
# [133] deldir_2.0-4              BiocParallel_1.36.0      
# [135] spatstat.geom_3.6-0       Matrix_1.6-5             
# [137] RcppHNSW_0.7.0            patchwork_1.3.2          
# [139] sparseMatrixStats_1.14.0  future_1.75.0            
# [141] shiny_1.11.1              ROCR_1.0-12              
# [143] igraph_2.1.4              memoise_2.0.1            