.libPaths()
#[1] "/micromamba/envs/seurat+milo.new/lib/R/library"

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

D20milo_k75<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_pca10_d10_prop0.1_2026_04_17.rds")
print("D20milo_k75@colData:")
head(D20milo_k75@colData)
da_results_D20<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_2026_04_17.rds")
print("da_results_D20:")
head(da_results_D20)
table(da_results_D20$celltypes_greek)
#prepare data for beeswarm plot
da_results_D20$"celltypes_greek"<-as.factor(da_results_D20$"celltypes_greek")
#check for NA values
na_values<-is.na(da_results_D20)
na_values<-which(na_values==TRUE)
#remove NA values
da_results_D20_filt<-da_results_D20[-na_values,]

# Make beeswarm graph
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_beeswarm_2026_04_17.pdf",width=15,height=10)
plotDAbeeswarm(da_results_D20_filt, group.by = "celltypes_greek")
dev.off()
da_results_D20_filt$celltypes_greek <- factor(da_results_D20_filt$celltypes_greek,
                              levels = rev(c("Monocytes","cDCs","pDCs","NK cells","CD2- γδ T cells","CD4+ αβ T cells","B cells","ASCs")))

CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_beeswarm_2026_04_17_v2.pdf", width=15, height=10)
plotDAbeeswarm(da_results_D20_filt, group.by = "celltypes_greek") +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
    axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
    axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
    axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
    legend.title = element_text(size = 20),                   # Customize legend title
    legend.text = element_text(size = 18)                    # Customize legend text
  ) 
dev.off()

CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf", width=15, height=10)
plotDAbeeswarm(da_results_D20_filt, group.by = "celltypes_greek") +
theme(
    axis.text.x = element_text(size = 24, color = "black", face = "bold"),  # Make x-axis text black
    axis.text.y = element_text(size = 24, color = "black", face = "bold"),  # Make y-axis text black
    axis.title.x = element_text(size = 22, color = "black", face = "bold"), # Make x-axis title black
    axis.title.y = element_text(size = 22, color = "black", face = "bold"), # Make y-axis title black
  )
  dev.off()
# Make beeswarm graph v2 with legend
#pdf(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_beeswarm_2026_04_17.pdf", width=15, height=10)
#plotDAbeeswarm(da_results_D20_filt, group.by = "celltypes_greek") +
  #theme(
   # axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
   # axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
   # axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
   # axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
   # legend.title = element_text(size = 20),                   # Customize legend title
   # legend.text = element_text(size = 18)                    # Customize legend text
  #) +
  #guides(color = guide_legend(title = "Cell Types"))         # Add legend for cell types
#dev.off()
#CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf", width=15, height=10)
#plotDAbeeswarm(da_results_D20_filt, group.by = "celltypes_greek") +
 # theme(
  #  axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
   #  axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
   # axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
    # axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
    # legend.title = element_text(size = 20),                   # Customize legend title
    # legend.text = element_text(size = 18)                    # Customize legend text
  #) +
  #guides(color = guide_legend(title = "Cell Types"))         # Add legend for cell types
#dev.off()
                  
#Let's furhter summarize only those neighborhoods with an annotated cell ID. Start by subsetting only non-Mixed neighborhoods, creating a column defining DA significance (FDR < 0.1), and creating a column indicating fold-change towards PP only or non-PP only enrichment:

da_D20_sum <- subset(da_results_D20_filt, celltypes_greek_fraction > 0.7)
da_D20_sum <- da_D20_sum %>% mutate(significant = case_when(SpatialFDR < 0.1 ~ 'Sig', SpatialFDR >= 0.1 ~ 'NonSig'))
da_D20_sum <- da_D20_sum %>% mutate(FC = case_when(logFC > 0 ~ '2DPI', logFC < 0 ~ '0DPI', logFC == 0 ~ 'Neutral'))
da_D20_sum$result <- paste(da_D20_sum$FC, da_D20_sum$significant, sep = '_')
da_D20_sum$result <- replace(da_D20_sum$result, da_D20_sum$result == '0DPI_NonSig', 'NonSig')
da_D20_sum$result <- replace(da_D20_sum$result, da_D20_sum$result == '2DPI_NonSig', 'NonSig')
# Reorder the levels of the 'result' factor
da_D20_sum$result <- factor(da_D20_sum$result, levels = c('0DPI_Sig', 'NonSig', '2DPI_Sig'))
da_D20_sum_table <- table(da_D20_sum$celltypes_greek,da_D20_sum$result)
write.table(da_D20_sum_table, 
            file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_sigNhoods_table_2026_04_17.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

percents <- prop.table(table(da_D20_sum$celltypes_greek,da_D20_sum$result),
                       margin = 1)
percents <- percents[,c('0DPI_Sig', 'NonSig', '2DPI_Sig')]
write.table(percents,file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_sigNhoods_percenttable_2026_04_17.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
percents <- t(percents) # transpose the table
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_stackedbarplot_2026_04_17.pdf", width=15, height=10)
# Adjust margins and mgp to move the x-axis title below the x-axis text
par(mar = c(10, 5, 4, 2) + 0.1, mgp = c(6, 1, 0))  # Increase bottom margin and adjust mgp
bp <- barplot(percents, # create stacked bar plot
              col = c('indianred4', 'grey60', 'darkslateblue'),
              xlab = "Cell type",  # X-axis title
              ylab = "Frequency of neighborhoods",  # Y-axis title
              border = NA,
              space = 0.05,
              names.arg = rep("", ncol(percents)))  # Suppress default x-axis labels

# Add custom x-axis labels at 45-degree angle
text(x = bp, y = -0.05, labels = colnames(percents), srt = 45, xpd = TRUE, adj = 1)
dev.off()

#make table of  # of neighbourhoods per celltypes with SpatialFDR < 0.1
sigNhoods<-table(da_results_D20$celltypes_greek[da_results_D20$SpatialFDR < 0.1])
sigNhoods<-as.data.frame(sigNhoods)
total_nhoods_per_celltype <- da_results_D20 %>%
  group_by(celltypes_greek) %>%
  summarise(total_Nhood = n())
colnames(sigNhoods) <- c("celltypes_greek", "Nhoods with SpatialFDR < 0.1")
sigNhoods <- merge(sigNhoods, total_nhoods_per_celltype, by="celltypes_greek")
colnames(sigNhoods) <- c("Cell type", "Nhoods with SpatialFDR < 0.1", "Total Nhoods per cell type")
sigNhoods$"sigNhoods/Total Nhoods" <- sigNhoods$"Nhoods with SpatialFDR < 0.1"/sigNhoods$"Total Nhoods per cell type"
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 0 per cell type
sig_nhoods_logFC_pos <- da_results_D20 %>%
  filter(SpatialFDR < 0.1 & logFC > 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 0")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 0 per cell type
sig_nhoods_logFC_neg <- da_results_D20 %>%
  filter(SpatialFDR < 0.1 & logFC < 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < 0")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_sigNhoods_celltype_2026_04_17.txt", sep="\t", row.names=FALSE)
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 1.5 per cell type
sig_nhoods_logFC_pos <- da_results_D20 %>%
  filter(SpatialFDR < 0.1 & logFC >= 1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 1.5")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 1.5 per cell type
sig_nhoods_logFC_neg <- da_results_D20 %>%
  filter(SpatialFDR < 0.1 & logFC <= -1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < 1.5")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
#replace NA values with 0
sigNhoods3[is.na(sigNhoods3)] <- 0
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D20_k75_DA_glmm_sigNhoods_celltype_2026_04_17_logFC15.txt", sep="\t", row.names=FALSE)




D80milo_k75<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_pca10_d10_prop0.1_2026_04_17.rds")
da_results_D80<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_2026_04_17.rds")
write.table(da_results_D80, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_2026_04_17.txt", sep="\t", row.names=FALSE)


#prepare data for beeswarm plot
da_results_D80$"celltypes_greek"<-as.factor(da_results_D80$"celltypes_greek")
#check for NA values
na_values<-is.na(da_results_D80)
na_values<-which(na_values==TRUE)
length(na_values)
#[1] 0


# Make beeswarm graph
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_beeswarm_2026_04_17.pdf", width=15, height=10)
plotDAbeeswarm(da_results_D80, group.by = "celltypes_greek")
dev.off()
# Make beeswarm graph v2
da_results_D80$celltypes_greek <- factor(da_results_D80$celltypes_greek,
                              levels = rev(c("Monocytes","cDCs","pDCs","NK cells","CD2- γδ T cells","CD4+ αβ T cells","B cells","ASCs")))
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_beeswarm_2026_04_17_v2.pdf", width=15, height=10)
plotDAbeeswarm(da_results_D80, group.by = "celltypes_greek") +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
    axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
    axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
    axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
    legend.title = element_text(size = 20),                   # Customize legend title
   legend.text = element_text(size = 18)                    # Customize legend text
  ) 
dev.off()

CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf",width=15,height=10)
plotDAbeeswarm(da_results_D80, group.by = "celltypes_greek")+
  theme(
    axis.text.x = element_text(size = 24, color = "black", face = "bold"),  # Make x-axis text black
    axis.text.y = element_text(size = 24, color = "black", face = "bold"),  # Make y-axis text black
    axis.title.x = element_text(size = 22, color = "black", face = "bold"), # Make x-axis title black
    axis.title.y = element_text(size = 22, color = "black", face = "bold"), # Make y-axis title black
  )
dev.off()
#CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf", width=15, height=10)
#plotDAbeeswarm(da_results_D80, group.by = "celltypes_greek") +
 # theme(
  #  axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
  #  axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
  # axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
  #  axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
  #  legend.title = element_text(size = 20),                   # Customize legend title
  #  legend.text = element_text(size = 18)                    # Customize legend text
  # ) +
 # guides(color = guide_legend(title = "Cell Types"))         # Add legend for cell types
#dev.off()

#Let's furhter summarize only those neighborhoods with an annotated cell ID. Start by subsetting only non-Mixed neighborhoods, creating a column defining DA significance (FDR < 0.1), and creating a column indicating fold-change towards PP only or non-PP only enrichment:

da_D80_sum <- subset(da_results_D80, celltypes_greek_fraction > 0.7)
da_D80_sum <- da_D80_sum %>% mutate(significant = case_when(SpatialFDR < 0.1 ~ 'Sig', SpatialFDR >= 0.1 ~ 'NonSig'))
da_D80_sum <- da_D80_sum %>% mutate(FC = case_when(logFC > 0 ~ '8DPI', logFC < 0 ~ '0DPI', logFC == 0 ~ 'Neutral'))
da_D80_sum$result <- paste(da_D80_sum$FC, da_D80_sum$significant, sep = '_')
da_D80_sum$result <- replace(da_D80_sum$result, da_D80_sum$result == '0DPI_NonSig', 'NonSig')
da_D80_sum$result <- replace(da_D80_sum$result, da_D80_sum$result == '8DPI_NonSig', 'NonSig')
# Reorder the levels of the 'result' factor
da_D80_sum$result <- factor(da_D80_sum$result, levels = c('0DPI_Sig', 'NonSig', '8DPI_Sig'))

da_D80_sum_table <- table(da_D80_sum$celltypes_greek,da_D80_sum$result)
write.table(da_D80_sum_table, 
            file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_sigNhoods_table_2026_04_17.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

percents <- prop.table(table(da_D80_sum$celltypes_greek,da_D80_sum$result),
                       margin = 1)
percents <- percents[,c('0DPI_Sig', 'NonSig', '8DPI_Sig')]
write.table(percents,file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_sigNhoods_percenttable_2026_04_17.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
percents <- t(percents) # transpose the table
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_stackedbarplot_2026_04_17.pdf", width=15, height=10)
# Adjust margins and mgp to move the x-axis title below the x-axis text
par(mar = c(10, 5, 4, 2) + 0.1, mgp = c(6, 1, 0))  # Increase bottom margin and adjust mgp
 bp <- barplot(percents, # create stacked bar plot
             col = c('indianred4', 'grey60', 'darkslateblue'),
              xlab = "Cell type",  # X-axis title
              ylab = "Frequency of neighborhoods",  # Y-axis title
              border = NA,
              space = 0.05,
              names.arg = rep("", ncol(percents)))  # Suppress default x-axis labels

# Add custom x-axis labels at 45-degree angle
text(x = bp, y = -0.05, labels = colnames(percents), srt = 45, xpd = TRUE, adj = 1)
dev.off()



#make table of  # of neighbourhoods per celltypes with SpatialFDR < 0.1
sigNhoods<-table(da_results_D80$celltypes_greek[da_results_D80$SpatialFDR < 0.1])
sigNhoods<-as.data.frame(sigNhoods)
# Calculate total Nhood per cell type
total_nhoods_per_celltype <- da_results_D80 %>%
  group_by(celltypes_greek) %>%
  summarise(total_Nhood = n())

colnames(sigNhoods) <- c("celltypes_greek", "Nhoods with SpatialFDR < 0.1")
sigNhoods <- merge(sigNhoods, total_nhoods_per_celltype, by="celltypes_greek")
colnames(sigNhoods) <- c("Cell type", "Nhoods with SpatialFDR < 0.1", "Total Nhoods per cell type")
sigNhoods$"sigNhoods/Total Nhoods" <- sigNhoods$"Nhoods with SpatialFDR < 0.1"/sigNhoods$"Total Nhoods per cell type"
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 0 per cell type
sig_nhoods_logFC_pos <- da_results_D80 %>%
  filter(SpatialFDR < 0.1 & logFC > 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 0")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 0 per cell type
sig_nhoods_logFC_neg <- da_results_D80 %>%
  filter(SpatialFDR < 0.1 & logFC < 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < 0")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_sigNhoods_celltype_2026_04_17.txt", sep="\t", row.names=FALSE)
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 1.5 per cell type
sig_nhoods_logFC_pos <- da_results_D80 %>%
  filter(SpatialFDR < 0.1 & logFC >= 1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 1.5")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 1.5 per cell type
sig_nhoods_logFC_neg <- da_results_D80 %>%
  filter(SpatialFDR < 0.1 & logFC <= -1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < 1.5")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
#replace NA values with 0
sigNhoods3[is.na(sigNhoods3)] <- 0
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D80_k75_DA_glmm_sigNhoods_celltype_2026_04_17_logFC15.txt", sep="\t", row.names=FALSE)



D82milo_k75<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_pca10_d10_prop0.1_2026_04_17.rds")
da_results_D82<-readRDS("/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_2026_04_17.rds")

#prepare data for beeswarm plot
da_results_D82$"celltypes_greek"<-as.factor(da_results_D82$"celltypes_greek")
#check for NA values
na_values<-is.na(da_results_D82)
na_values<-which(na_values==TRUE)
length(na_values)
#[1] 1688
#remove NA values
da_results_D82_filt<-da_results_D82[-na_values,]

# Make beeswarm graph
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_beeswarm_2026_04_17.pdf",width=15,height=10)
plotDAbeeswarm(da_results_D82_filt, group.by = "celltypes_greek")
dev.off()
# Make beeswarm graph v2
da_results_D82_filt$celltypes_greek <- factor(da_results_D82_filt$celltypes_greek,
                              levels = rev(c("Monocytes","cDCs","pDCs","NK cells","CD2- γδ T cells","CD4+ αβ T cells","B cells","ASCs")))
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_beeswarm_2026_04_17_v2.pdf", width=15, height=10)
plotDAbeeswarm(da_results_D82_filt, group.by = "celltypes_greek") +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
    axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
    axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
    axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
    legend.title = element_text(size = 20),                   # Customize legend title
    legend.text = element_text(size = 18)                    # Customize legend text
  )
dev.off()


CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf",width=15,height=10)
plotDAbeeswarm(da_results_D82_filt, group.by = "celltypes_greek")+
  theme(
    axis.text.x = element_text(size = 24, color = "black", face = "bold"),  # Make x-axis text black
    axis.text.y = element_text(size = 24, color = "black", face = "bold"),  # Make y-axis text black
     axis.title.x = element_text(size = 22, color = "black", face = "bold"), # Make x-axis title black
     axis.title.y = element_text(size = 22, color = "black", face = "bold"), # Make y-axis title black
  )
dev.off()

#CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_beeswarm_2026_04_17_v3.pdf", width=15, height=10)
#plotDAbeeswarm(da_results_D82_filt, group.by = "celltypes_greek") +
 # theme(
 #   axis.text.x = element_text(size = 20, color = "black"),  # Make x-axis text black
  #  axis.text.y = element_text(size = 20, color = "black"),  # Make y-axis text black
  #  axis.title.x = element_text(size = 22, color = "black"), # Make x-axis title black
  #  axis.title.y = element_text(size = 22, color = "black"), # Make y-axis title black
  #  legend.title = element_text(size = 20),                   # Customize legend title
  #  legend.text = element_text(size = 18)                    # Customize legend text
  #) +
 # guides(color = guide_legend(title = "Cell Types"))         # Add legend for cell types
#dev.off()

#Let's furhter summarize only those neighborhoods with an annotated cell ID. Start by subsetting only non-Mixed neighborhoods, creating a column defining DA significance (FDR < 0.1), and creating a column indicating fold-change towards PP only or non-PP only enrichment:

da_D82_sum <- subset(da_results_D82_filt, celltypes_greek_fraction > 0.7)
da_D82_sum <- da_D82_sum %>% mutate(significant = case_when(SpatialFDR < 0.1 ~ 'Sig', SpatialFDR >= 0.1 ~ 'NonSig'))
da_D82_sum <- da_D82_sum %>% mutate(FC = case_when(logFC > 0 ~ '8DPI', logFC < 0 ~ '2DPI', logFC == 0 ~ 'Neutral'))
da_D82_sum$result <- paste(da_D82_sum$FC, da_D82_sum$significant, sep = '_')
da_D82_sum$result <- replace(da_D82_sum$result, da_D82_sum$result == '2DPI_NonSig', 'NonSig')
da_D82_sum$result <- replace(da_D82_sum$result, da_D82_sum$result == '8DPI_NonSig', 'NonSig')
# Reorder the levels of the 'result' factor
da_D82_sum$result <- factor(da_D82_sum$result, levels = c('2DPI_Sig', 'NonSig', '8DPI_Sig'))

da_D82_sum_table <- table(da_D82_sum$celltypes_greek,da_D82_sum$result)
write.table(da_D82_sum_table, 
            file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_sigNhoods_table_2026_04_17.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

percents <- prop.table(table(da_D82_sum$celltypes_greek,da_D82_sum$result),
                       margin = 1)
percents <- percents[,c('2DPI_Sig', 'NonSig', '8DPI_Sig')]
write.table(percents,file = "/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_sigNhoods_percenttable_2026_04_17.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
percents <- t(percents) # transpose the table
CairoPDF(file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_stackedbarplot_2026_04_17.pdf", width=15, height=10)
# Adjust margins and mgp to move the x-axis title below the x-axis text
par(mar = c(10, 5, 4, 2) + 0.1, mgp = c(6, 1, 0))  # Increase bottom margin and adjust mgp
bp <- barplot(percents, # create stacked bar plot
              col = c('indianred4', 'grey60', 'darkslateblue'),
              xlab = "Cell type",  # X-axis title
              ylab = "Frequency of neighborhoods",  # Y-axis title
              border = NA,
              space = 0.05,
              names.arg = rep("", ncol(percents)))  # Suppress default x-axis labels

# Add custom x-axis labels at 45-degree angle
text(x = bp, y = -0.05, labels = colnames(percents), srt = 45, xpd = TRUE, adj = 1)
dev.off()

#make table of  # of neighbourhoods per celltypes with SpatialFDR < 0.1
sigNhoods<-table(da_results_D82$celltypes_greek[da_results_D82$SpatialFDR < 0.1])
sigNhoods<-as.data.frame(sigNhoods)
# Calculate total Nhood per cell type
total_nhoods_per_celltype <- da_results_D82 %>%
  group_by(celltypes_greek) %>%
  summarise(total_Nhood = n())

colnames(sigNhoods) <- c("celltypes_greek", "Nhoods with SpatialFDR < 0.1")
sigNhoods <- merge(sigNhoods, total_nhoods_per_celltype, by="celltypes_greek")
colnames(sigNhoods) <- c("Cell type", "Nhoods with SpatialFDR < 0.1", "Total Nhoods per cell type")
sigNhoods$"sigNhoods/Total Nhoods" <- sigNhoods$"Nhoods with SpatialFDR < 0.1"/sigNhoods$"Total Nhoods per cell type"
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 0 per cell type
sig_nhoods_logFC_pos <- da_results_D82 %>%
  filter(SpatialFDR < 0.1 & logFC > 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 0")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 0 per cell type
sig_nhoods_logFC_neg <- da_results_D82 %>%
  filter(SpatialFDR < 0.1 & logFC < 0) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < 0")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_sigNhoods_celltype_2026_04_17.txt", sep="\t", row.names=FALSE)
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC > 1.5 per cell type
sig_nhoods_logFC_pos <- da_results_D82 %>%
  filter(SpatialFDR < 0.1 & logFC >= 1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_pos = n())
colnames(sig_nhoods_logFC_pos) <- c("Cell type", "sigNhoods with logFC > 1.5")
# Calculate number of Nhoods with SpatialFDR < 0.1 and logFC < 1.5 per cell type
sig_nhoods_logFC_neg <- da_results_D82 %>%
  filter(SpatialFDR < 0.1 & logFC <= -1.5) %>%
  group_by(celltypes_greek) %>%
  summarise(Nhoods_with_SpatialFDR_and_logFC_neg = n())
colnames(sig_nhoods_logFC_neg) <- c("Cell type", "sigNhoods with logFC < -1.5")
sigNhoods2 <-left_join(sigNhoods,sig_nhoods_logFC_neg, by="Cell type")
sigNhoods3 <- left_join(sigNhoods2,sig_nhoods_logFC_pos, by="Cell type")
#replace NA values with 0
sigNhoods3[is.na(sigNhoods3)] <- 0
write.table(sigNhoods3, file="/scRNAseq/Sal_5pigs_2026/DA/milo_D82_k75_DA_glmm_sigNhoods_celltype_2026_04_17_LogFC15.txt", sep="\t", row.names=FALSE)

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
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] Cairo_1.6-2   cowplot_1.1.3 ggplot2_3.5.1 dplyr_1.1.4   miloR_2.1.3
# [6] edgeR_4.0.16  limma_3.58.1

# loaded via a namespace (and not attached):
#  [1] tidyselect_1.2.1            viridisLite_0.4.2
#  [3] vipor_0.4.7                 farver_2.1.2
#  [5] viridis_0.6.5               bitops_1.0-9
#  [7] ggraph_2.2.1                fastmap_1.2.0
#  [9] SingleCellExperiment_1.24.0 RCurl_1.98-1.16
# [11] pracma_2.4.4                tweenr_2.0.3
# [13] rsvd_1.0.5                  lifecycle_1.0.4
# [15] statmod_1.5.0               magrittr_2.0.3
# [17] compiler_4.3.3              rlang_1.1.4
# [19] tools_4.3.3                 igraph_2.1.2
# [21] labeling_0.4.3              S4Arrays_1.2.1
# [23] graphlayouts_1.2.1          DelayedArray_0.28.0
# [25] RColorBrewer_1.1-3          abind_1.4-8
# [27] BiocParallel_1.36.0         numDeriv_2016.8-1.1
# [29] withr_3.0.2                 purrr_1.0.2
# [31] BiocGenerics_0.48.1         grid_4.3.3
# [33] polyclip_1.10-7             stats4_4.3.3
# [35] beachmat_2.18.1             colorspace_2.1-1
# [37] scales_1.3.0                gtools_3.9.5
# [39] MASS_7.3-60                 SummarizedExperiment_1.32.0
# [41] cli_3.6.3                   crayon_1.5.3
# [43] generics_0.1.3              ggbeeswarm_0.7.2
# [45] cachem_1.1.0                ggforce_0.4.2
# [47] stringr_1.5.1               zlibbioc_1.48.2
# [49] parallel_4.3.3              XVector_0.42.0
# [51] matrixStats_1.4.1           vctrs_0.6.5
# [53] Matrix_1.6-5                BiocSingular_1.18.0
# [55] IRanges_2.36.0              patchwork_1.3.0
# [57] S4Vectors_0.40.2            BiocNeighbors_1.20.2
# [59] ggrepel_0.9.6               irlba_2.3.5.1
# [61] beeswarm_0.4.0              locfit_1.5-9.10
# [63] tidyr_1.3.1                 glue_1.8.0
# [65] codetools_0.2-20            stringi_1.8.4
# [67] gtable_0.3.6                GenomeInfoDb_1.38.8
# [69] GenomicRanges_1.54.1        ScaledMatrix_1.10.0
# [71] munsell_0.5.1               tibble_3.2.1
# [73] pillar_1.10.0               GenomeInfoDbData_1.2.11
# [75] R6_2.5.1                    tidygraph_1.3.1
# [77] lattice_0.22-6              Biobase_2.62.0
# [79] memoise_2.0.1               Rcpp_1.0.13-1
# [81] gridExtra_2.3               SparseArray_1.2.4
# [83] MatrixGenerics_1.14.0       pkgconfig_2.0.3