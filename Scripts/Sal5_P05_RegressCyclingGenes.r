.libPaths()
#[1] "/micromamba/envs/seurat+milo/lib/R/library"


{set.seed(123)
library(Seurat)
library(SeuratExtendData)
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
library(findPC)
library(tidyr)
library(cowplot)
library(scCustomize)
library(SingleCellExperiment)}



# Load in SCE object from Scran env
All_sce.2<-readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_ScranNorm_SCE_2026_04_15.rds")
pca_dims <- reducedDim(All_sce.2, "PCA")
# Make seurat object with normalized counts
All<- as.Seurat(All_sce.2, counts="counts", data = "logcounts")
# rename the PCA dimensions to "pca" so seurat recognizes them
All[["pca"]] <- CreateDimReducObject(embeddings= pca_dims, key = "PC_", assay = DefaultAssay(All))
All[["PCA"]] <- NULL
All$time_point <- sub('.*_(D\\d+)_.*', '\\1', All$SampleID)
#Change assay name from originalexp to RNA
All[['RNA']] = All[['originalexp']]
DefaultAssay(All) <- "RNA"
All[['originalexp']] = NULL
#add scale.data to seurat object
All<- ScaleData(All)
#Add Variable Features
# Assuming nfeatures is a vector of feature names
nfeatures=dget("/scRNAseq/Sal_5pigs_2026/Normalization/HVGs_scran.R")
VariableFeatures(All) <- nfeatures
#Save seurat object with doublet score in metadata
saveRDS(All,file ="/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_ScranNorm_seurat_2026_04_15.rds")
 
All<-readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_ScranNorm_seurat_2026_04_15.rds")
#Make new seurat object with normalized counts ONLY
# get normalized counts
cDat <- as.matrix(GetAssayData(object = All, layer = 'data'))
#create phenotype data
pDat <-data.frame(barcode = colnames(All))
pDat$SampleID <- sapply(strsplit(pDat$barcode, "_"), function(x) paste(x[1:4], collapse = "_"))
pDat$BarBak <- pDat$barcode
pDat <- pDat %>% separate(BarBak, c("Sample","Animal"))
#Warning message:Expected 2 pieces. Additional pieces discarded in 142577 rows [1, 2, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
rownames(pDat) <- pDat$barcode
pDat$GenesDetected <- colSums(cDat!=0)
pDat$UmiSums<- colSums(cDat)
All_scran <- CreateSeuratObject(counts = cDat, meta.data = pDat) # create Seurat object of counts & pheno data

Idents(object = All_scran)<-"SampleID"


pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scranNorm_postnorm_scaledplot.pdf",width=10, height= 10)
scatter_plot <- QC_Plot_UMIvsGene(seurat_object = All_scran)
# Modify the x-axis limits using ggplot2
scatter_plot <- scatter_plot + xlim(0, 40000) +ylim(0, 6000) + theme(aspect.ratio = 1)
scatter_plot
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scranNorm_postnorm_plot.pdf",width=10, height= 10)
QC_Plot_UMIvsGene(seurat_object = All_scran)
dev.off()

annotKey<-read.csv("/Annotation_files/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.csv",header=T)
table(annotKey$gene_biotype)
#select genes gene_biotype = "ribozyme"
riboGenes<-annotKey[annotKey$gene_biotype == "ribozyme",]
ribo_gene_list <- riboGenes$gene_name
mitoGenes<-read.csv("/Annotation_files/mitogenes.csv",header=T)
mito_gene_list <- mitoGenes$gene
All_scran<- Add_Mito_Ribo(object = All_scran, species = "other", mito_features = mito_gene_list, ribo_features = ribo_gene_list)

pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/scranNorm_PostNorm_gene_UMI_mito.pdf", width=25, height= 10)
QC_Plots_Combined_Vln(All_scran, group.by = "SampleID",plot_boxplot = TRUE,x_lab_rotate = TRUE)
dev.off()

pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/PostScran_genes_histo.pdf",width=25, height= 10)
QC_Histogram(seurat_object = All, features = "nFeature_RNA",split.by = "SampleID",high_cutoff = 400)
dev.off()

# expression distribution after Scran normalization sample 10,000 genes
ScranNorm_geneExp = as.vector(All_scran[['RNA']]$counts) %>% sample(10000)
# remove non-expressed genes
ScranNorm_geneExp = ScranNorm_geneExp[ScranNorm_geneExp != 0]
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/postScranNorm_gene_exp_histo.pdf")
hist(ScranNorm_geneExp)
dev.off()



#try with findPC
    # Get the standard deviation of each PC
    stdev <- All[["pca"]]@stdev
    # Calculate the variance explained by each PC
    var_explained <- stdev^2 / sum(stdev^2) * 100

# pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/Sal5_P04_SCRAN_findPC.pdf", width = 30, height = 10)
# findPC(sdev = stdev,number = c(10:20),method = 'all',figure = T)
# dev.off()

# pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/Sal5_P04_SCRAN_findPCs6t20.pdf", width = 30, height = 10)
# findPC(sdev = stdev,number = c(6:20),method = 'all',figure = T)
# dev.off()

#Regressed out cell cycle genes & redo PCA
All<-readRDS("/scRNAseq/Sal_5pigs_2026/Orig_files/5pigs_salmonella_Filtered_ScranNorm_seurat_2026_04_15.rds")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# get cell cycle scores
All <- CellCycleScoring(All, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regressed out cycling genes
All  <- ScaleData(All , vars.to.Regressed = c("S.Score", "G2M.Score"), features = rownames(All))
All<- RunPCA(All, npcs = 100, verbose = TRUE)
saveRDS(All, file ="/scRNAseq/Sal_5pigs_2026/Normalization/5pigs_salmonella_Filtered_ScranNorm_Regressed_PCA_2026_04_15.rds")
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/ScranNorm_Regressed_PCA_elbow.pdf")
ElbowPlot(All, ndims = 40)
dev.off()
    stdev <- All[["pca"]]@stdev
    # Calculate the variance explained by each PC
    var_explained <- stdev^2 / sum(stdev^2) * 100
# Find the PC that explains less than 1% of the variance
pc_below_1_percent <- which(var_explained < 1)[1]  # Get the first PC below #Only keep the first 30 PCs 
var_explained <- var_explained[1:30]

    # Plot the percentage of variance explained with line at 1% of variance explained with legend
    pdf("/scRNAseq/Sal_5pigs_2026/Normalization/ScranNorm_Regressed_PCA_ElbowPlot_variance_explained.pdf")
    plot(var_explained, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained", xaxt = "n")
    axis(1, at = 1:length(var_explained), labels = 1:length(var_explained))
    abline(h = 1, col = "red")
    abline(v = pc_below_1_percent, col = "red")
    legend("topright", legend = c("1% of variance explained"), col = c("red"), lty = 1)
    dev.off()    
# Create a data frame with the variance explained and the principal component number
pca_pct_variance <- data.frame(PC = 1:length(var_explained), variance = var_explained)
# Use PCAtools to find the elbow point
print("Use PCAtools to find the elbow point")
chosen_elbow <- PCAtools::findElbowPoint(pca_pct_variance$variance)
# Plot PCAtools elbow
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/ScranNorm_Regressed_PCA_PCAtools_elbow.pdf")
pca_pct_variance %>% 
  ggplot(aes(PC, variance)) + 
  geom_point() + 
  geom_vline(xintercept = chosen_elbow) + 
  ylab("% of variance explained")
  dev.off()
# How many PCs explain more than 1% of the variance in the entire dataset.
print("table(pca_pct_variance$variance > 1)")
table(pca_pct_variance$variance > 1)
##FALSE  TRUE
##  86     14
#Which PCs explain less than 1% of the variance in the entire dataset.
print("table(pca_pct_variance$variance < 1)")
table(pca_pct_variance$variance < 1)
##FALSE  TRUE
##  14     86
#PC14 explains 1.0378312% of the variance in the entire dataset & has Preceding residual=2 .
#try with findPC
pdf(file="/scRNAseq/Sal_5pigs_2026/Normalization/ScranNorm_Regressed_PCA_findPC.pdf", width = 30, height = 10)
findPC(sdev = stdev,number = c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),method = 'all',figure = T)
dev.off()

#Want low number for Preceding residual & second derivative
findPCs<-findPC(sdev = stdev,number = c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),method = 'all',figure = F)
class(findPCs)
write.table(findPCs,file="/scRNAseq/Sal_5pigs_2026/Normalization/ScranNorm_Regressed_PCA_findPCs.txt",sep="\t",quote=F,row.names=T)
print("findPC(sdev = stdev,number = c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),method = 'all',aggregate = 'median',figure = F)")
findPC(sdev = stdev,number = c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),method = 'all',aggregate = 'median',figure = F)

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
#  [3] tidyr_1.3.1                 findPC_1.0                 
#  [5] scater_1.30.1               scuttle_1.12.0             
#  [7] DESeq2_1.42.1               harmony_1.2.3              
#  [9] Rcpp_1.0.14                 writexl_1.5.1              
# [11] openxlsx_4.2.8              DropletUtils_1.22.0        
# [13] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [15] Biobase_2.62.0              GenomicRanges_1.54.1       
# [17] GenomeInfoDb_1.38.8         IRanges_2.36.0             
# [19] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# [21] MatrixGenerics_1.14.0       matrixStats_1.5.0          
# [23] patchwork_1.3.0             magrittr_2.0.3             
# [25] dplyr_1.1.4                 ggplot2_3.5.2              
# [27] Matrix_1.6-5                SeuratExtendData_0.2.0     
# [29] Seurat_5.2.1                SeuratObject_5.0.2         
# [31] sp_2.2-0                   

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.3.3            
#   [3] later_1.4.2               prismatic_1.1.2          
#   [5] bitops_1.0-9              tibble_3.2.1             
#   [7] R.oo_1.27.0               polyclip_1.10-7          
#   [9] janitor_2.2.1             fastDummies_1.7.5        
#  [11] lifecycle_1.0.4           edgeR_4.0.16             
#  [13] globals_0.16.3            lattice_0.22-6           
#  [15] MASS_7.3-60               limma_3.58.1             
#  [17] plotly_4.10.4             httpuv_1.6.15            
#  [19] sctransform_0.4.1         spam_2.11-1              
#  [21] zip_2.3.2                 spatstat.sparse_3.1-0    
#  [23] reticulate_1.41.0.1       pbapply_1.7-2            
#  [25] RColorBrewer_1.1-3        lubridate_1.9.4          
#  [27] abind_1.4-8               zlibbioc_1.48.2          
#  [29] Rtsne_0.17                purrr_1.0.4              
#  [31] R.utils_2.13.0            RCurl_1.98-1.17          
#  [33] circlize_0.4.16           GenomeInfoDbData_1.2.11  
#  [35] ggrepel_0.9.6             irlba_2.3.5.1            
#  [37] listenv_0.9.1             spatstat.utils_3.1-2     
#  [39] goftest_1.2-3             RSpectra_0.16-2          
#  [41] spatstat.random_3.3-2     dqrng_0.4.1              
#  [43] fitdistrplus_1.2-2        parallelly_1.42.0        
#  [45] DelayedMatrixStats_1.24.0 codetools_0.2-20         
#  [47] DelayedArray_0.28.0       shape_1.4.6.1            
#  [49] tidyselect_1.2.1          farver_2.1.2             
#  [51] viridis_0.6.5             ScaledMatrix_1.10.0      
#  [53] spatstat.explore_3.3-4    jsonlite_2.0.0           
#  [55] BiocNeighbors_1.20.2      progressr_0.15.1         
#  [57] ggridges_0.5.6            survival_3.8-3           
#  [59] tools_4.3.3               ica_1.0-3                
#  [61] glue_1.8.0                gridExtra_2.3            
#  [63] SparseArray_1.2.4         HDF5Array_1.30.1         
#  [65] withr_3.0.2               fastmap_1.2.0            
#  [67] rhdf5filters_1.14.1       rsvd_1.0.5               
#  [69] digest_0.6.37             timechange_0.3.0         
#  [71] R6_2.6.1                  mime_0.13                
#  [73] ggprism_1.0.5             colorspace_2.1-1         
#  [75] scattermore_1.2           tensor_1.5               
#  [77] spatstat.data_3.1-4       R.methodsS3_1.8.2        
#  [79] generics_0.1.3            data.table_1.17.0        
#  [81] httr_1.4.7                htmlwidgets_1.6.4        
#  [83] S4Arrays_1.2.1            uwot_0.2.3               
#  [85] pkgconfig_2.0.3           gtable_0.3.6             
#  [87] lmtest_0.9-40             XVector_0.42.0           
#  [89] htmltools_0.5.8.1         dotCall64_1.2            
#  [91] scales_1.4.0              png_0.1-8                
#  [93] snakecase_0.11.1          spatstat.univar_3.1-2    
#  [95] reshape2_1.4.4            nlme_3.1-167             
#  [97] GlobalOptions_0.1.2       zoo_1.8-13               
#  [99] rhdf5_2.46.1              stringr_1.5.1            
# [101] KernSmooth_2.23-26        vipor_0.4.7              
# [103] parallel_4.3.3            miniUI_0.1.1.1           
# [105] ggrastr_1.0.2             pillar_1.10.2            
# [107] grid_4.3.3                vctrs_0.6.5              
# [109] RANN_2.6.2                promises_1.3.2           
# [111] BiocSingular_1.18.0       beachmat_2.18.1          
# [113] xtable_1.8-4              cluster_2.1.8.1          
# [115] paletteer_1.6.0           beeswarm_0.4.0           
# [117] cli_3.6.5                 locfit_1.5-9.12          
# [119] compiler_4.3.3            rlang_1.1.6              
# [121] crayon_1.5.3              future.apply_1.11.3      
# [123] labeling_0.4.3            rematch2_2.1.2           
# [125] forcats_1.0.0             ggbeeswarm_0.7.2         
# [127] plyr_1.8.9                stringi_1.8.7            
# [129] viridisLite_0.4.2         deldir_2.0-4             
# [131] BiocParallel_1.36.0       lazyeval_0.2.2           
# [133] spatstat.geom_3.3-5       PCAtools_2.14.0          
# [135] RcppHNSW_0.6.0            sparseMatrixStats_1.14.0 
# [137] future_1.34.0             Rhdf5lib_1.24.2          
# [139] statmod_1.5.0             shiny_1.10.0             
# [141] ROCR_1.0-11               igraph_2.1.4             