#======================================#
#         LOAD LIBRARIES
#======================================#

.libPaths()
#[1] "/micromamba/envs/seurat+milo/lib/R/library"
{set.seed(123)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(scCustomize)}



# Load in SCE object from scDoublet env
All_sce.2<-readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_QC_IDed_Doublets_sce_2026_04_15.rds")
# Make seurat object with doublet score in metadata
All_seurat <- CreateSeuratObject(counts = assays(All_sce.2)$counts, meta.data = as.data.frame(colData(All_sce.2)))
#23155 features across 182215 sampless
#Save seurat object with doublet score in metadata
saveRDS(All_seurat,file ="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_QC_IDed_Doublets_seurat_2026_04_15.rds")

#filter by scDblFinder.class(several metrics are used to give cell a "class")
Idents(object = All_seurat)<-"scDblFinder.class"
All<-subset(x = All_seurat, idents ="singlet")
# 23155 features across 151808 samples


# Make scCustomize QC plots
Idents(object = All)<-"SampleID"

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/gene_UMI_mito_doublets_removed.pdf", width=25, height= 10)
QC_Plots_Combined_Vln(All, group.by = "SampleID",plot_boxplot = TRUE,x_lab_rotate = TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/Complexity_violin_doublets_removed.pdf", width =15, height= 15)
QC_Plots_Complexity(seurat_object = All, high_cutoff = 0.8,pt.size = 0, y_axis_log = TRUE, plot_median = TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/mito_histo_doublets_removed.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "percent_mito_ribo",split.by = "SampleID",high_cutoff = 15)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/genes_histo_doublets_removed.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "nFeature_RNA",split.by = "SampleID",high_cutoff = 400)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/Umis_histo_doublets_removed.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "nCount_RNA",split.by = "SampleID",high_cutoff = 300)
dev.off()

{raw_median_stats <- Median_Stats(seurat_object = All, group_by_var = "SampleID")
raw_median_stats<-as.data.frame(raw_median_stats)
raw_meta<-All@meta.data
ncells<-table(All$SampleID)
ncells<-as.data.frame(ncells)
colnames(ncells)[1]<- "SampleID"
colnames(ncells)[2]<- "nCells"
raw_all<- dplyr::left_join(raw_median_stats, ncells, by = "SampleID")
head(raw_all)
raw_all_s<-write.csv(raw_all,file="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_PreProcessing_QC_stats_2026_04_15.csv")
raw_all_s2<-write.csv(raw_all,file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/5pigs_PreProcessing_QC_stats_doublets_removed_2026_04_15.csv")}

#Add info to  metadata to use as filters
{All$PassViability=All$percent_mito_ribo < 15
All$PassGenesDet=All$nFeature_RNA> 400
All$PassLibSize=All$nCount_RNA > 300
All$PassBarcodeFreq=All$DuplicatedBarcodes==FALSE
All$PassAll= All$PassViability & All$PassGenesDet & All$PassLibSize & All$PassBarcodeFreq}





Idents(object = All)<-"PassAll"
All<-subset(x = All, idents ="TRUE")
# 23155 features across 84569 samples

#Check that all is filtered
met_df<-All@meta.data
dplyr::filter(met_df, PassAll== "FLASE")
#<0 rows> (or 0-length row.names)


# Now save filtered seurat object
saveRDS(All,file ="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_QC_IDed_Doublets_Filtered_seurat_2026_04_15.rds")
{All_counts<-All[["RNA"]]$counts
All_metadata<-All@meta.data
All_sce<-SingleCellExperiment(assays=list(counts=All_counts), colData=All_metadata)
sce_s<- saveRDS(All_sce,"/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_SCE_2026_04_15.rds")}


pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/scatter_Prenorm_plot.pdf",width=30, height= 10)
P1<-QC_Plot_UMIvsGene(seurat_object = All, low_cutoff_gene = 400,high_cutoff_gene = 5500, low_cutoff_UMI = 300)
P2<-QC_Plot_UMIvsGene(seurat_object = All, meta_gradient_name = "percent_mito", low_cutoff_gene = 600,high_cutoff_gene = 5500, high_cutoff_UMI = 45000)
P1|P2
dev.off()

{raw_median_stats <- Median_Stats(seurat_object = All, group_by_var = "SampleID")
raw_median_stats<-as.data.frame(raw_median_stats)
raw_meta<-All@meta.data
ncells<-table(All$SampleID)
ncells<-as.data.frame(ncells)
colnames(ncells)[1]<- "SampleID"
colnames(ncells)[2]<- "nCells"
raw_all<- dplyr::left_join(raw_median_stats, ncells, by = "SampleID")
head(raw_all)
raw_all_s<-write.csv(raw_all,file="/scRNAseq/Sal_5pigs_2026/Orig_files/PostProcessing_QC_stats.csv")
raw_all_s2<-write.csv(raw_all,file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/PostProcessing_QC_stats.csv")}

#Plot original expression distribution
raw_geneExp = as.vector(All[['RNA']]$counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/Raw_gene_exp_histo.pdf")
hist(raw_geneExp)
dev.off()

Idents(object = All)<-"SampleID"

pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/gene_UMI_mito.pdf", width=25, height= 10)
QC_Plots_Combined_Vln(All, group.by = "SampleID",plot_boxplot = TRUE,x_lab_rotate = TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/Complexity_violin.pdf", width =15, height= 15)
QC_Plots_Complexity(seurat_object = All, high_cutoff = 0.8,pt.size = 0, y_axis_log = TRUE, plot_median = TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/mito_histo.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "percent_mito_ribo",split.by = "SampleID",high_cutoff = 15)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/genes_histo.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "nFeature_RNA",split.by = "SampleID",high_cutoff = 400)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/PostProcessing_Plots/Umis_histo.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "nCount_RNA",split.by = "SampleID",high_cutoff = 300)
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#  [1] scCustomize_3.0.1           SingleCellExperiment_1.24.0
#  [3] SummarizedExperiment_1.32.0 Biobase_2.62.0             
#  [5] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
#  [7] IRanges_2.36.0              S4Vectors_0.40.2           
#  [9] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
# [11] matrixStats_1.5.0           dplyr_1.1.4                
# [13] ggplot2_3.5.2               Matrix_1.6-5               
# [15] Seurat_5.2.1                SeuratObject_5.0.2         
# [17] sp_2.2-0                   

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      shape_1.4.6.1           jsonlite_2.0.0         
#   [4] magrittr_2.0.3          ggbeeswarm_0.7.2        spatstat.utils_3.1-2   
#   [7] farver_2.1.2            GlobalOptions_0.1.2     zlibbioc_1.48.2        
#  [10] vctrs_0.6.5             ROCR_1.0-11             spatstat.explore_3.3-4 
#  [13] paletteer_1.6.0         RCurl_1.98-1.17         janitor_2.2.1          
#  [16] forcats_1.0.0           htmltools_0.5.8.1       S4Arrays_1.2.1         
#  [19] SparseArray_1.2.4       sctransform_0.4.1       parallelly_1.42.0      
#  [22] KernSmooth_2.23-26      htmlwidgets_1.6.4       ica_1.0-3              
#  [25] plyr_1.8.9              lubridate_1.9.4         plotly_4.10.4          
#  [28] zoo_1.8-13              igraph_2.1.4            mime_0.13              
#  [31] lifecycle_1.0.4         pkgconfig_2.0.3         R6_2.6.1               
#  [34] fastmap_1.2.0           snakecase_0.11.1        GenomeInfoDbData_1.2.11
#  [37] fitdistrplus_1.2-2      future_1.34.0           shiny_1.10.0           
#  [40] digest_0.6.37           colorspace_2.1-1        rematch2_2.1.2         
#  [43] patchwork_1.3.0         tensor_1.5              prismatic_1.1.2        
#  [46] RSpectra_0.16-2         irlba_2.3.5.1           labeling_0.4.3         
#  [49] progressr_0.15.1        timechange_0.3.0        spatstat.sparse_3.1-0  
#  [52] httr_1.4.7              polyclip_1.10-7         abind_1.4-8            
#  [55] compiler_4.3.3          withr_3.0.2             fastDummies_1.7.5      
#  [58] MASS_7.3-60             DelayedArray_0.28.0     tools_4.3.3            
#  [61] vipor_0.4.7             lmtest_0.9-40           beeswarm_0.4.0         
#  [64] httpuv_1.6.15           future.apply_1.11.3     goftest_1.2-3          
#  [67] glue_1.8.0              nlme_3.1-167            promises_1.3.2         
#  [70] grid_4.3.3              Rtsne_0.17              cluster_2.1.8.1        
#  [73] reshape2_1.4.4          generics_0.1.3          gtable_0.3.6           
#  [76] spatstat.data_3.1-4     tidyr_1.3.1             data.table_1.17.0      
#  [79] XVector_0.42.0          spatstat.geom_3.3-5     RcppAnnoy_0.0.22       
#  [82] ggrepel_0.9.6           RANN_2.6.2              pillar_1.10.2          
#  [85] stringr_1.5.1           ggprism_1.0.5           spam_2.11-1            
#  [88] RcppHNSW_0.6.0          later_1.4.2             circlize_0.4.16        
#  [91] splines_4.3.3           lattice_0.22-6          survival_3.8-3         
#  [94] deldir_2.0-4            tidyselect_1.2.1        miniUI_0.1.1.1         
#  [97] pbapply_1.7-2           gridExtra_2.3           scattermore_1.2        
# [100] stringi_1.8.7           lazyeval_0.2.2          codetools_0.2-20       
# [103] tibble_3.2.1            cli_3.6.5               uwot_0.2.3             
# [106] xtable_1.8-4            reticulate_1.41.0.1     Rcpp_1.0.14            
# [109] globals_0.16.3          spatstat.random_3.3-2   png_0.1-8              
# [112] ggrastr_1.0.2           spatstat.univar_3.1-2   parallel_4.3.3         
# [115] dotCall64_1.2           bitops_1.0-9            listenv_0.9.1          
# [118] viridisLite_0.4.2       scales_1.4.0            ggridges_0.5.6         
# [121] purrr_1.0.4             crayon_1.5.3            rlang_1.1.6            
# [124] cowplot_1.1.3          