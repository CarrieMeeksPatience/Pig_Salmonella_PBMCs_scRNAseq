#======================================#
#         LOAD LIBRARIES
#======================================#

.libPaths()
#[1] "/micromamba/envs/seurat+milo.new/lib/R/library"
{library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(DropletUtils)
library(openxlsx)
library(writexl)
library(harmony)
library(DESeq2)
library(scater)
library(tidyr)
library(cowplot)
library(scCustomize)}

#======================================#
#              QC 
#======================================#
setwd('/PIPseq_PIPseeker/ALL/results')
data_dir.2 <-c(
  S01_852_D0_S_1 = '/PIPseq_PIPseeker/ALL/results/852_0D_S_1_merged/filtered_matrix/sensitivity_4',
  S02_852_D2_S_2 = '/PIPseq_PIPseeker/ALL/results/852_2D_S_2/filtered_matrix/sensitivity_3',
  S03_852_D8_S_2 = '/PIPseq_PIPseeker/ALL/results/852_8D_S_2/filtered_matrix/sensitivity_3',
  S04_854_D0_S_2 = "/PIPseq_PIPseeker/ALL/results/854_0D_S_2_Novogene/filtered_matrix/sensitivity_5",
  S05_854_D2_S_2 = "/PIPseq_PIPseeker/ALL/results/854_2D_S_2_merged/filtered_matrix/sensitivity_5",
  S06_854_D8_S_2 = "/PIPseq_PIPseeker/ALL/results/854_8D_S_2_merged/filtered_matrix/sensitivity_5",
  S07_864_D0_S_1 = "/PIPseq_PIPseeker/ALL/results/864_0D_S_1_merged/filtered_matrix/sensitivity_4",
  S08_864_D2_S_2 = "/PIPseq_PIPseeker/ALL/results/864_2D_S_2_merged/filtered_matrix/sensitivity_3",
  S09_864_D8_S_2 = "/PIPseq_PIPseeker/ALL/results/864_8D_S_2_merged/filtered_matrix/sensitivity_3",
  S10_842_D0_S_1 = "/PIPseq_PIPseeker/ALL/results/842_0D_S_1-results/filtered_matrix/sensitivity_5",
  S11_842_D2_S_1 = "/PIPseq_PIPseeker/ALL/results/842_2D_S_1_merged-results/filtered_matrix/sensitivity_5",
  S12_842_D8_S_2 = "/PIPseq_PIPseeker/ALL/results/842_8D_S_2-results/filtered_matrix/sensitivity_5",
  S13_853_D0_S_1 = "/PIPseq_PIPseeker/ALL/results/853_0D_S_1_merged2-results/filtered_matrix/sensitivity_5",
  S14_853_D2_S_1 = "/PIPseq_PIPseeker/ALL/results/853_2D_S_1-results/filtered_matrix/sensitivity_5",
  S15_853_D8_S_2 = "/PIPseq_PIPseeker/ALL/results/853_8D_S_2_merged-results/filtered_matrix/sensitivity_5"
)

#read data directory
library_id<- c("S01_852_D0_S_1","S02_852_D2_S_2","S03_852_D8_S_2","S04_854_D0_S_2","S05_854_D2_S_2","S06_854_D8_S_2","S07_864_D0_S_1","S08_864_D2_S_2","S09_864_D8_S_2","S10_842_D0_S_1","S11_842_D2_S_1","S12_842_D8_S_2","S13_853_D0_S_1","S14_853_D2_S_1","S15_853_D8_S_2")
scRNA_data <- Read10X(data.dir = data_dir.2)
seurat_object = CreateSeuratObject(counts = scRNA_data)
#create count matrix
cDat <- as.matrix(GetAssayData(object = seurat_object, layer = 'counts'))
dim(cDat) 
#[1]  25880 146293
#remove non-expressed genes and more than 0 transcript across all cells
keep <- rowSums(cDat) > 0
cDat <- cDat[keep,]
dim(cDat) 
#[1]   21486 146293

# create feature data
fDat <- data.frame(ID = rownames(cDat))
rownames(fDat) <- fDat$ID
con <- gzfile(file.path(data_dir.2[1], "features.tsv.gz"))
ssc_genes <- read.delim(con, sep = "\t", header = FALSE, as.is = TRUE)[, 1:2]
colnames(ssc_genes) <- c("PIP", "Dashed")
ssc_genes$Dashed<- gsub("_", "-", ssc_genes$Dashed, perl = TRUE); head(ssc_genes, n = 3)
#                 PIP  Dashed
# 1 ENSSSCG00000037372     TBP
# 2 ENSSSCG00000027257   PSMB1
# 3 ENSSSCG00000029697 FAM120B
ssc_genes$Symbol <- sub("_.*", "", ssc_genes$PIP)
ssc_genes$EnsemblID <- sub(".*_", "", ssc_genes$PIP)
ssc_genes$Duplicated <- duplicated(ssc_genes$Symbol) | duplicated(ssc_genes$Symbol, fromLast = TRUE)
ssc_genes$Name <- ifelse(ssc_genes$Duplicated == "TRUE", ssc_genes$EnsemblID, ssc_genes$Symbol)
fDat <- merge(fDat, ssc_genes, by.x ="row.names", by.y = "Dashed", all.x =TRUE, sort = FALSE)
rownames(fDat) <- fDat[, 1]; fDat <- fDat[, -1]

#Read in jayne's mitochrondia genes
mitoGenes<-read.csv("/Annotation_files/mitogenes.csv",header=T)

fDat$Mitochondrial <- fDat$EnsemblID %in% mitoGenes$x
table(fDat$Mitochondrial)
#FALSE  TRUE
#21449     37


#create phenotype data
pDat <-data.frame(barcode = colnames(seurat_object))
pDat$SampleID <- sapply(strsplit(pDat$barcode, "_"), function(x) paste(x[1:4], collapse = "_"))
pDat$BarBak <- pDat$barcode
pDat <- pDat %>% separate(BarBak, c("Sample","Animal"))
#Warning message:Expected 2 pieces. Additional pieces discarded in 142577 rows [1, 2, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
rownames(pDat) <- pDat$barcode
##Check for duplicated barcodes(UMIs)
pDat$DuplicatedBarcodes <- duplicated(rownames(pDat)) | duplicated(rownames(pDat), fromLast = TRUE)
table(pDat$DuplicatedBarcodes)
#True  FALSE
#0     146293 

#Make seurat object to create QC plots
All <- CreateSeuratObject(counts = cDat, meta.data = pDat)
All
#21486 features across 146293 samples

##Make scCustomize QC plots
annotKey<-read.csv("/Annotation_files/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.csv",header=T)
table(annotKey$gene_biotype)
#select genes gene_biotype = "ribozyme"
riboGenes<-annotKey[annotKey$gene_biotype == "ribozyme",]
ribo_gene_list <- riboGenes$gene_id
mito_gene_list <- mitoGenes$gene
All<- Add_Mito_Ribo(object = All, species = "other", mito_features = mito_gene_list, ribo_features = ribo_gene_list)
All <- Add_Cell_Complexity(object = All)


#Make QC stats table
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
raw_all_s2<-write.csv(raw_all,file="/scRNAseq/Sal_5pigs_2026/PreProcessing_Plots/5pigs_PreProcessing_QC_stats_2026_04_15.csv")}


#Save Seurat object
All_s <- saveRDS(All,"/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_seurat_QC_2026_04_15.rds")

#Make SCE object too!
# Get counts data
{All_counts<-All[["RNA"]]$counts
All_metadata<-All@meta.data
All_sce<-SingleCellExperiment(assays=list(counts=All_counts), colData=All_metadata)
All_sce_s<-saveRDS(All_sce,"/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_SCE_QC_2026_04_15.rds")}

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
#  [1] scCustomize_3.0.1           cowplot_1.1.3              
#  [3] tidyr_1.3.1                 scater_1.30.1              
#  [5] scuttle_1.12.0              DESeq2_1.42.1              
#  [7] harmony_1.2.3               Rcpp_1.0.14                
#  [9] writexl_1.5.1               openxlsx_4.2.8             
# [11] DropletUtils_1.22.0         SingleCellExperiment_1.24.0
# [13] SummarizedExperiment_1.32.0 Biobase_2.62.0             
# [15] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
# [17] IRanges_2.36.0              S4Vectors_0.40.2           
# [19] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
# [21] matrixStats_1.5.0           patchwork_1.3.0            
# [23] magrittr_2.0.3              dplyr_1.1.4                
# [25] ggplot2_3.5.2               Matrix_1.6-5               
# [27] Seurat_5.2.1                SeuratObject_5.0.2         
# [29] sp_2.2-0                   

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.3.3            
#   [3] later_1.4.2               bitops_1.0-9             
#   [5] tibble_3.2.1              R.oo_1.27.0              
#   [7] polyclip_1.10-7           janitor_2.2.1            
#   [9] fastDummies_1.7.5         lifecycle_1.0.4          
#  [11] edgeR_4.0.16              globals_0.16.3           
#  [13] lattice_0.22-6            MASS_7.3-60              
#  [15] limma_3.58.1              plotly_4.10.4            
#  [17] httpuv_1.6.15             sctransform_0.4.1        
#  [19] spam_2.11-1               zip_2.3.2                
#  [21] spatstat.sparse_3.1-0     reticulate_1.41.0.1      
#  [23] pbapply_1.7-2             RColorBrewer_1.1-3       
#  [25] lubridate_1.9.4           abind_1.4-8              
#  [27] zlibbioc_1.48.2           Rtsne_0.17               
#  [29] purrr_1.0.4               R.utils_2.13.0           
#  [31] RCurl_1.98-1.17           circlize_0.4.16          
#  [33] GenomeInfoDbData_1.2.11   ggrepel_0.9.6            
#  [35] irlba_2.3.5.1             listenv_0.9.1            
#  [37] spatstat.utils_3.1-2      goftest_1.2-3            
#  [39] RSpectra_0.16-2           spatstat.random_3.3-2    
#  [41] dqrng_0.4.1               fitdistrplus_1.2-2       
#  [43] parallelly_1.42.0         DelayedMatrixStats_1.24.0
#  [45] codetools_0.2-20          DelayedArray_0.28.0      
#  [47] shape_1.4.6.1             tidyselect_1.2.1         
#  [49] farver_2.1.2              viridis_0.6.5            
#  [51] ScaledMatrix_1.10.0       spatstat.explore_3.3-4   
#  [53] jsonlite_2.0.0            BiocNeighbors_1.20.2     
#  [55] progressr_0.15.1          ggridges_0.5.6           
#  [57] survival_3.8-3            tools_4.3.3              
#  [59] ica_1.0-3                 glue_1.8.0               
#  [61] gridExtra_2.3             SparseArray_1.2.4        
#  [63] HDF5Array_1.30.1          withr_3.0.2              
#  [65] fastmap_1.2.0             rhdf5filters_1.14.1      
#  [67] digest_0.6.37             rsvd_1.0.5               
#  [69] timechange_0.3.0          R6_2.6.1                 
#  [71] mime_0.13                 ggprism_1.0.5            
#  [73] colorspace_2.1-1          scattermore_1.2          
#  [75] tensor_1.5                spatstat.data_3.1-4      
#  [77] R.methodsS3_1.8.2         generics_0.1.3           
#  [79] data.table_1.17.0         httr_1.4.7               
#  [81] htmlwidgets_1.6.4         S4Arrays_1.2.1           
#  [83] uwot_0.2.3                pkgconfig_2.0.3          
#  [85] gtable_0.3.6              lmtest_0.9-40            
#  [87] XVector_0.42.0            htmltools_0.5.8.1        
#  [89] dotCall64_1.2             scales_1.4.0             
#  [91] png_0.1-8                 snakecase_0.11.1         
#  [93] spatstat.univar_3.1-2     reshape2_1.4.4           
#  [95] nlme_3.1-167              GlobalOptions_0.1.2      
#  [97] zoo_1.8-13                rhdf5_2.46.1             
#  [99] stringr_1.5.1             KernSmooth_2.23-26       
# [101] vipor_0.4.7               parallel_4.3.3           
# [103] miniUI_0.1.1.1            ggrastr_1.0.2            
# [105] pillar_1.10.2             grid_4.3.3               
# [107] vctrs_0.6.5               RANN_2.6.2               
# [109] promises_1.3.2            BiocSingular_1.18.0      
# [111] beachmat_2.18.1           xtable_1.8-4             
# [113] cluster_2.1.8.1           paletteer_1.6.0          
# [115] beeswarm_0.4.0            cli_3.6.5                
# [117] locfit_1.5-9.12           compiler_4.3.3           
# [119] rlang_1.1.6               crayon_1.5.3             
# [121] future.apply_1.11.3       rematch2_2.1.2           
# [123] forcats_1.0.0             ggbeeswarm_0.7.2         
# [125] plyr_1.8.9                stringi_1.8.7            
# [127] viridisLite_0.4.2         deldir_2.0-4             
# [129] BiocParallel_1.36.0       lazyeval_0.2.2           
# [131] spatstat.geom_3.3-5       RcppHNSW_0.6.0           
# [133] sparseMatrixStats_1.14.0  future_1.34.0            
# [135] Rhdf5lib_1.24.2           statmod_1.5.0            
# [137] shiny_1.10.0              ROCR_1.0-11              
# [139] igraph_2.1.4             