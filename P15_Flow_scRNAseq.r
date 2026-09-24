.libPaths()
# [1] "/micromamba/envs/Emmeans/lib/R/library"
{library(lme4)
library(emmeans)
library(stats)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(hrbrthemes)
library(gridExtra)
library(rstatix)
library(ggpubr)
}
#datasets

Flow <- read.csv("/scRNAseq/Sal_5pigs_2026/scSal_Flow/FlowProportons_CellConcentration_2025_02_13_scSalPaper_V2_Rready.csv", header = T)
#Make column "DPI" from "Sample"
Flow <- Flow %>%
  mutate(DPI = gsub(".*_(D\\d+)_.*", "\\1", Sample))
table(Flow$DPI)
Flow$DPI <- gsub("D0", "0 DPI", Flow$DPI)
Flow$DPI <- gsub("D2", "2 DPI", Flow$DPI)
Flow$DPI <- gsub("D8", "8 DPI", Flow$DPI)
# Set DPI as a factor
Flow$DPI <- as.factor(Flow$DPI)

#Make column "Animal" from "Sample", the 1st 3 digits of Sample are the animal number
Flow <- Flow %>%
  mutate(Animal = gsub("^(\\d{3}).*", "\\1", Sample))
#Make column "celltype" from "cell.type"
table(Flow$Animal)
# Set Animal as a factor
Flow$Animal <- as.factor(Flow$Animal)
head(Flow)
#Replace Flow$celltype==ASCs(Dump) with Flow$celltype==ASCs
Flow$celltype <- gsub("ASCs\\(Dump\\)", "ASCs", Flow$celltype)
table(Flow$celltype)
#For  Flow$celltype=="", replace with "CD3e- CD172a+ CD8a+"
Flow$celltype <- gsub("^$", "CD3e- CD172a+ CD8a+", Flow$celltype)
table(Flow$celltype)
Flow$celltype <- gsub("CD2\\- GD T-cells", "CD2\\- gd T cells", Flow$celltype)
Flow$celltype <- gsub("CD4\\+ AB T-cells", "CD4\\+ ab T cells", Flow$celltype)
Flow$celltype <- gsub("monocytes CD4\\-", "Monocytes", Flow$celltype)
Flow$celltype <- gsub("pDCs & cDCs CD4\\+", "pDCs & cDCs", Flow$celltype)
Flow$celltype <- gsub("B-cells", "B cells", Flow$celltype)

#Change Flow$Percent_of_PBMCs name to Flow$Portion_of_PBMCs
Flow <- Flow %>%
  rename(Portion_of_PBMCs = Percent_of_PBMCs)
#make histogram of celltype_per_ul to see if it is right-skewed for each cell type
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_histogram.pdf", width = 12, height = 8)
ggplot(Flow, aes(x = celltype_per_ul)) +
  geom_histogram(bins = 8, fill = "blue", alpha = 0.7) +
  facet_wrap(~ celltype, scales = "free") +
  theme_minimal() +
  labs(
    title = "Cells per microliter; Cell-type distribution across 15 samples from Flow cytometry",
    x = "Cells per microliter",
    y = "Number of samples"
  )
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_boxplot.pdf", width = 12, height = 8)
ggplot(Flow, aes(x = DPI, y = celltype_per_ul)) +
  geom_boxplot(bins = 8, fill = "blue", alpha = 0.7) +     geom_jitter(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Cells per microliter; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Cells per microliter"
  )
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_line.pdf", width = 12, height = 8)
ggplot(Flow, aes(x = DPI, y = celltype_per_ul, group= Animal)) +
  geom_line( alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Cells per microliter; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Cells per microliter"
  )
dev.off()
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_line_V2.pdf", width = 12, height = 8)
ggplot(Flow, aes(x = DPI, y = celltype_per_ul, group = Animal,
    color = Animal)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Cells per microliter; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Cells per microliter"
  )
dev.off()

Flow_V2 <- Flow %>%
  filter( !celltype %in% c("CD3e- CD172a+ CD8a+", "PBMCs") )
cell_type <- unique(Flow_V2$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit <- list()
model_results <- list()
model_normality <- list()
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_residuals_D0_D2_D8.pdf",
  width = 12,
  height = 8
)
par(mfrow = c(1, 2))
for (cell_type in unique(Flow_V2$celltype)) {
  cell_type_data <- Flow_V2 %>%
    filter(celltype == cell_type)

  fit <- lmer(
    celltype_per_ul ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit[[cell_type]] <- fit
  model_normality[[cell_type]] <- shapiro.test(residuals(fit))
  model_results[[cell_type]] <- pairs(emmeans(fit, ~ DPI), adjust = "holm")
  model_fit[[cell_type]] <- fit
  model_results[[cell_type]] <-
    pairs(emmeans(fit, ~ DPI), adjust = "holm")

  hist(
    residuals(fit),
    main = paste("Residuals:", cell_type),
    xlab = "Residuals"
  )

  car::qqPlot(
    residuals(fit),
    main = paste("Normal Q-Q plot:", cell_type),
    id = FALSE
  )
}

dev.off()
capture.output(
  lapply(model_normality, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_shapiroTest_summaries_D0_D2_D8.txt"
)

# All Flow cell types show non-normality of residuals, so we will use non-parametric tests for pairwise comparisons between DPI groups.
Flow_D02 <- Flow %>%
  filter(
    DPI %in% c("0 DPI", "2 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+", "PBMCs")
  )
cell_type <- unique(Flow_D02$celltype)

model_fit_D02 <- list()
model_results_D02 <- list()

for (cell_type in unique(Flow_D02$celltype)) {
  cell_type_data <- Flow_D02 %>%
    filter(celltype == cell_type)

  D02_fit <- lmer(
    celltype_per_ul ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit_D02[[cell_type]] <- D02_fit

  emm <- emmeans(D02_fit, ~ DPI)

  model_results_D02[[cell_type]] <- contrast(
    emm,
    method = list("D2-D0" = c(-1, 1)),
    adjust = "holm"
  )
  
}

capture.output(
  lapply(model_fit_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D0_D2.txt"
)

capture.output(
  lapply(model_results_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D2.txt"
)
saveRDS(model_fit_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D0_D2.rds")
saveRDS(model_results_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")

Flow_D28 <- Flow %>%
  filter(
    DPI %in% c("2 DPI", "8 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+", "PBMCs")
  )
cell_type <- unique(Flow_D28$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D28 <- list()
model_results_D28 <- list()


for (cell_type in unique(Flow_D28$celltype)) {
  cell_type_data <- Flow_D28 %>%
    filter(celltype == cell_type)

  D28_fit <- lmer(
    celltype_per_ul ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D28[[cell_type]] <- D28_fit
  emm <- emmeans(D28_fit, ~ DPI)
  model_results_D28[[cell_type]] <- contrast(
    emm,
    method = list("D8-D2" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D2_D8.txt"
)

capture.output(
  lapply(model_results_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D2_D8.txt"
)
saveRDS(model_fit_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D2_D8.rds")
saveRDS(model_results_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")

Flow_D08 <- Flow %>%
  filter(
    DPI %in% c("0 DPI", "8 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+", "PBMCs")
  )
cell_type <- unique(Flow_D08$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D08 <- list()
model_results_D08 <- list()

for (cell_type in unique(Flow_D08$celltype)) {
  cell_type_data <- Flow_D08 %>%
    filter(celltype == cell_type)

  D08_fit <- lmer(
    celltype_per_ul ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D08[[cell_type]] <- D08_fit
  emm <- emmeans(D08_fit, ~ DPI)
  model_results_D08[[cell_type]] <- contrast(
    emm,
    method = list("D8-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D0_D8.txt"
)

capture.output(
  lapply(model_results_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D8.txt"
)
saveRDS(model_fit_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_summaries_D0_D8.rds")
saveRDS(model_results_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
##Now repeat with "Portion_of_PBMCs"  replacing "celltype_per_ul" in the above code with "Portion_of_PBMCs" and save the results to a new file
Flow_PofPBMCS <- Flow %>%
  filter( !celltype %in% c("Neutrophils", "PBMCs") )
#make histogram of Portion_of_PBMCs to see if it is right-skewed for each cell type
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_histogram.pdf", width = 12, height = 8)
ggplot(Flow_PofPBMCS, aes(x = Portion_of_PBMCs)) +
  geom_histogram(bins = 8, fill = "blue", alpha = 0.7) +
  facet_wrap(~ celltype, scales = "free") +
  theme_minimal() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from Flow cytometry",
    x = "Portion of PBMCs",
    y = "Number of samples"
  )
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_boxplot.pdf", width = 12, height = 8)
ggplot(Flow_PofPBMCS, aes(x = DPI, y = Portion_of_PBMCs)) +
  geom_boxplot(bins = 8, fill = "blue", alpha = 0.7) +     geom_jitter(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_line.pdf", width = 12, height = 8)
ggplot(Flow_PofPBMCS, aes(x = DPI, y = Portion_of_PBMCs, group= Animal)) +
  geom_line( alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_line_V2.pdf", width = 12, height = 8)
ggplot(Flow_PofPBMCS, aes(x = DPI, y = Portion_of_PBMCs, group = Animal,
    color = Animal)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from Flow cytometry",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()

rm(Flow_V2,model_fit,model_results,model_normality,fit)
Flow_V2 <- Flow %>%
  filter( !celltype %in% c("CD3e- CD172a+ CD8a+", "Neutrophils", "PBMCs") )
cell_type <- unique(Flow_V2$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit <- list()
model_results <- list()
model_normality <- list()
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_residuals_D0_D2_D8.pdf",
  width = 12,
  height = 8
)
par(mfrow = c(1, 2))
for (cell_type in unique(Flow_V2$celltype)) {
  cell_type_data <- Flow_V2 %>%
    filter(celltype == cell_type)

  fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit[[cell_type]] <- fit
  model_normality[[cell_type]] <- shapiro.test(residuals(fit))
  model_results[[cell_type]] <- pairs(emmeans(fit, ~ DPI), adjust = "holm")
  model_fit[[cell_type]] <- fit
  model_results[[cell_type]] <-
    pairs(emmeans(fit, ~ DPI), adjust = "holm")

  hist(
    residuals(fit),
    main = paste("Residuals:", cell_type),
    xlab = "Residuals"
  )

  car::qqPlot(
    residuals(fit),
    main = paste("Normal Q-Q plot:", cell_type),
    id = FALSE
  )
}

dev.off()
capture.output(
  lapply(model_normality, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_shapiroTest_summaries_D0_D2_D8.txt"
)

# All Flow cell types show non-normality of residuals, so we will use non-parametric tests for pairwise comparisons between DPI groups.
rm(Flow_D02,model_fit_D02,model_results_D02,D02_fit)
Flow_D02 <- Flow %>%
  filter(
    DPI %in% c("0 DPI", "2 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+","Neutrophils", "PBMCs")
  )
cell_type <- unique(Flow_D02$celltype)

model_fit_D02 <- list()
model_results_D02 <- list()

for (cell_type in unique(Flow_D02$celltype)) {
  cell_type_data <- Flow_D02 %>%
    filter(celltype == cell_type)

  D02_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit_D02[[cell_type]] <- D02_fit
  
  emm <- emmeans(D02_fit, ~ DPI)

  model_results_D02[[cell_type]] <- contrast(
    emm,
    method = list("D2-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D2.txt"
)

capture.output(
  lapply(model_results_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.txt"
)

saveRDS(model_fit_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D2.rds")
saveRDS(model_results_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")
rm(Flow_D28,model_fit_D28,model_results_D28,D28_fit)
Flow_D28 <- Flow %>%
  filter(
    DPI %in% c("2 DPI", "8 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+","Neutrophils", "PBMCs")
  )
cell_type <- unique(Flow_D28$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D28 <- list()
model_results_D28 <- list()


for (cell_type in unique(Flow_D28$celltype)) {
  cell_type_data <- Flow_D28 %>%
    filter(celltype == cell_type)

  D28_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D28[[cell_type]] <- D28_fit
  emm <- emmeans(D28_fit, ~ DPI)
  model_results_D28[[cell_type]] <- contrast(
    emm,
    method = list("D8-D2" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D2_D8.txt"
)

capture.output(
  lapply(model_results_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.txt"
)
saveRDS(model_fit_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D2_D8.rds")
saveRDS(model_results_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")

rm(Flow_D08,model_fit_D08,model_results_D08,D08_fit)
Flow_D08 <- Flow %>%
  filter(
    DPI %in% c("0 DPI", "8 DPI"),
    !celltype %in% c("CD3e- CD172a+ CD8a+","Neutrophils", "PBMCs")
  )
cell_type <- unique(Flow_D08$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D08 <- list()
model_results_D08 <- list()

for (cell_type in unique(Flow_D08$celltype)) {
  cell_type_data <- Flow_D08 %>%
    filter(celltype == cell_type)

  D08_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D08[[cell_type]] <- D08_fit
  emm <- emmeans(D08_fit, ~ DPI)
  model_results_D08[[cell_type]] <- contrast(
    emm,
    method = list("D8-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D8.txt"
)

capture.output(
  lapply(model_results_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.txt"
)

saveRDS(model_fit_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D8.rds")
saveRDS(model_results_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")

## Now we will anaylze the scRNAseq data 
scRNAseq_nCells <- read.table("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Animal_Number_Prop_long.txt", sep = "\t", header = T)
scRNAseq_nCells$Sample_Celltype <- paste(scRNAseq_nCells$Sample, scRNAseq_nCells$Celltype, sep = "_")
head(scRNAseq_nCells)
scRNAseq_prop <- read.table("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Sample100_Prop_long.txt", sep = "\t", header = T)
scRNAseq_prop$Sample_Celltype <- paste(scRNAseq_prop$Sample, scRNAseq_prop$Celltype, sep = "_")
head(scRNAseq_prop)
#merge the two dataframes by Sample_Celltype
scRNAseq <- merge(scRNAseq_nCells, scRNAseq_prop, by = "Sample_Celltype")
head(scRNAseq)
scRNAseq <- scRNAseq %>%
  rename(celltype = Celltype.x)
scRNAseq <- scRNAseq %>%
  rename(Sample = Sample.x)
scRNAseq$Celltype.y <- NULL
scRNAseq$Sample.y <- NULL
scRNAseq$Technology <- "scRNAseq"
#scRNAseq$Sample values are formated as "###_D#".Pull the "D#" part to create the DPI column.
scRNAseq$DPI <- sub(".*_(D\\d+)$", "\\1", scRNAseq$Sample)
table(scRNAseq$DPI)
scRNAseq$DPI <- gsub("D0", "0 DPI", scRNAseq$DPI)
scRNAseq$DPI <- gsub("D2", "2 DPI", scRNAseq$DPI)
scRNAseq$DPI <- gsub("D8", "8 DPI", scRNAseq$DPI)
scRNAseq$DPI <- as.factor(scRNAseq$DPI)
#scRNAseq$Sample values are formated as "###_D#".Pull the "###" part to create the Animal column.
scRNAseq$Animal <- gsub("^(\\d+)_D.*", "\\1", scRNAseq$Sample)
table(scRNAseq$Animal)
#Change scRNAseq$Proportion.of.Sample name to scRNAseq$Portion_of_PBMCs
scRNAseq <- scRNAseq %>%
  rename(Portion_of_PBMCs = Proportion.of.Sample)
head(scRNAseq)
table(scRNAseq$celltype)
scRNAseq_pDC_cDC <- scRNAseq %>%
  filter(celltype %in% c("pDCs", "cDCs")) %>%
  group_by(Sample, DPI, Animal) %>%
  summarise(
    celltype = "pDCs & cDCs",
    nCells = sum(nCells, na.rm = TRUE),
    Portion_of_PBMCs = sum(Portion_of_PBMCs, na.rm = TRUE),
    .groups = "drop"
  )
print("head of scRNAseq_pDC_cDC")
head(scRNAseq_pDC_cDC)
#Add scRNAseq_pDC_cDC to scRNAseq
scRNAseq <- bind_rows(scRNAseq, scRNAseq_pDC_cDC)
#make histogram of Portion_of_PBMCs to see if it is right-skewed for each cell type
scRNAseq_V2 <- scRNAseq %>%
  filter( !celltype %in% c("Unknown", "Total Cells in Sample") )
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_histogram.pdf", width = 12, height = 8)
ggplot(scRNAseq_V2, aes(x = Portion_of_PBMCs)) +
  geom_histogram(bins = 8, fill = "blue", alpha = 0.7) +
  facet_wrap(~ celltype, scales = "free") +
  theme_minimal() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from scRNAseq",
    x = "Portion of PBMCs",
    y = "Number of samples"
  )
dev.off()


pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_boxplot.pdf", width = 12, height = 8)
ggplot(scRNAseq_V2, aes(x = DPI, y = Portion_of_PBMCs)) +
  geom_boxplot(bins = 8, fill = "blue", alpha = 0.7) +     geom_jitter(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from scRNAseq",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_line.pdf", width = 12, height = 8)
ggplot(scRNAseq_V2, aes(x = DPI, y = Portion_of_PBMCs, group= Animal)) +
  geom_line( alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from scRNAseq",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_line_V2.pdf", width = 12, height = 8)
ggplot(scRNAseq_V2, aes(x = DPI, y = Portion_of_PBMCs, group = Animal,
    color = Animal)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from scRNAseq",
    x = "DPI",
    y = "Portion of PBMCs"
  )
dev.off()

#First run analysis on scRNAseq$celltype=="Total Cells in Sample" only, to see if there is a significant difference in the total number of cells between DPI groups. 
scRNAseq_V3 <- scRNAseq %>%
  filter( celltype %in% c("Total Cells in Sample") )
cell_type <- unique(scRNAseq_V3$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit <- list()
model_results <- list()
model_normality <- list()
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_residuals_D0_D2_D8.pdf",
  width = 12,
  height = 8
)
par(mfrow = c(1, 2))
for (cell_type in unique(scRNAseq_V3$celltype)) {
  cell_type_data <- scRNAseq_V3 %>%
    filter(celltype == cell_type)

  fit <- lmer(
    nCells ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit[[cell_type]] <- fit
  model_normality[[cell_type]] <- shapiro.test(residuals(fit))
  model_results[[cell_type]] <- pairs(emmeans(fit, ~ DPI), adjust = "holm")
  model_fit[[cell_type]] <- fit
  model_results[[cell_type]] <-
    pairs(emmeans(fit, ~ DPI), adjust = "holm")

  hist(
    residuals(fit),
    main = paste("Residuals:", cell_type),
    xlab = "Residuals"
  )

  car::qqPlot(
    residuals(fit),
    main = paste("Normal Q-Q plot:", cell_type),
    id = FALSE
  )
}

dev.off()
capture.output(
  lapply(model_normality, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_shapiroTest_summaries_D0_D2_D8.txt"
)

# All Flow cell types show non-normality of residuals, so we will use non-parametric tests for pairwise comparisons between DPI groups.
scRNAseq_V3_D02 <- scRNAseq %>%
  filter(
    DPI %in% c("0 DPI", "2 DPI"),
    celltype %in% c("Total Cells in Sample")
  )
cell_type <- unique(scRNAseq_V3_D02$celltype)

model_fit_D02 <- list()
model_results_D02 <- list()

for (cell_type in unique(scRNAseq_V3_D02$celltype)) {
  cell_type_data <- scRNAseq_V3_D02 %>%
    filter(celltype == cell_type)

  D02_fit <- lmer(
    nCells ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit_D02[[cell_type]] <- D02_fit

  model_results_D02[[cell_type]] <- pairs(emmeans(D02_fit, ~ DPI), adjust = "holm")
  emm <- emmeans(D02_fit, ~ DPI)
  model_results_D02[[cell_type]] <- contrast(
    emm,
    method = list("D2-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D0_D2.txt"
)

capture.output(
  lapply(model_results_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D0_D2.txt"
)

saveRDS(model_fit_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D0_D2.rds")
saveRDS(model_results_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")

scRNAseq_V3_D28 <- scRNAseq %>%
  filter(
    DPI %in% c("2 DPI", "8 DPI"),
    celltype %in% c("Total Cells in Sample")
  )
cell_type <- unique(scRNAseq_V3_D28$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D28 <- list()
model_results_D28 <- list()



for (cell_type in unique(scRNAseq_V3_D28$celltype)) {
  cell_type_data <- scRNAseq_V3_D28 %>%
    filter(celltype == cell_type)

  D28_fit <- lmer(
    nCells ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D28[[cell_type]] <- D28_fit
  emm <- emmeans(D28_fit, ~ DPI)
  model_results_D28[[cell_type]] <- contrast(
    emm,
    method = list("D8-D2" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D2_D8.txt"
)

capture.output(
  lapply(model_results_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D2_D8.txt"
)

saveRDS(model_fit_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D2_D8.rds")
saveRDS(model_results_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")

scRNAseq_V3_D08 <- scRNAseq %>%
  filter(
    DPI %in% c("0 DPI", "8 DPI"),
    celltype %in% c("Total Cells in Sample")
  )
cell_type <- unique(scRNAseq_V3_D08$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D08 <- list()
model_results_D08 <- list()

for (cell_type in unique(scRNAseq_V3_D08$celltype)) {
  cell_type_data <- scRNAseq_V3_D08 %>%
    filter(celltype == cell_type)

  D08_fit <- lmer(
    nCells ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  emm <- emmeans(D08_fit, ~ DPI)
  model_results_D08[[cell_type]] <- contrast(
    emm,
    method = list("D8-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D0_D8.txt"
)

capture.output(
  lapply(model_results_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D0_D8.txt"
)

saveRDS(model_fit_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_summaries_D0_D8.rds")
saveRDS(model_results_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_nCells_per_Sample_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
##Now repeat with "Portion_of_PBMCs"  replacing "celltype_per_ul" in the above code with "Portion_of_PBMCs" and save the results to a new file
scRNAseq_V2 <- scRNAseq %>%
  filter( !celltype %in% c("Unknown", "Total Cells in Sample") )
cell_type <- unique(scRNAseq_V2$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit <- list()
model_results <- list()
model_normality <- list()
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_residuals_D0_D2_D8.pdf",
  width = 12,
  height = 8
)
par(mfrow = c(1, 2))
for (cell_type in unique(scRNAseq_V2$celltype)) {
  cell_type_data <- scRNAseq_V2 %>%
    filter(celltype == cell_type)

  fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit[[cell_type]] <- fit
  model_normality[[cell_type]] <- shapiro.test(residuals(fit))
  model_results[[cell_type]] <- pairs(emmeans(fit, ~ DPI), adjust = "holm")
  model_fit[[cell_type]] <- fit
  model_results[[cell_type]] <-
    pairs(emmeans(fit, ~ DPI), adjust = "holm")

  hist(
    residuals(fit),
    main = paste("Residuals:", cell_type),
    xlab = "Residuals"
  )

  car::qqPlot(
    residuals(fit),
    main = paste("Normal Q-Q plot:", cell_type),
    id = FALSE
  )
}

dev.off()
capture.output(
  lapply(model_normality, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_shapiroTest_summaries_D0_D2_D8.txt"
)

# All Flow cell types show non-normality of residuals, so we will use non-parametric tests for pairwise comparisons between DPI groups.
scRNAseq_D02 <- scRNAseq %>%
  filter(
    DPI %in% c("0 DPI", "2 DPI"),
    !celltype %in% c("Unknown", "Total Cells in Sample") 
  )
cell_type <- unique(scRNAseq_D02$celltype)

model_fit_D02 <- list()
model_results_D02 <- list()

for (cell_type in unique(scRNAseq_D02$celltype)) {
  cell_type_data <- scRNAseq_D02 %>%
    filter(celltype == cell_type)

  D02_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )
  model_fit_D02[[cell_type]] <- D02_fit
  model_results_D02[[cell_type]] <- pairs(emmeans(D02_fit, ~ DPI), adjust = "holm")
  emm <- emmeans(D02_fit, ~ DPI)
  model_results_D02[[cell_type]] <- contrast(
    emm,
    method = list("D2-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D2.txt"
)

capture.output(
  lapply(model_results_D02, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.txt"
)

saveRDS(model_fit_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D2.rds")
saveRDS(model_results_D02, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")
scRNAseq_D28 <- scRNAseq %>%
  filter(
    DPI %in% c("2 DPI", "8 DPI"),
    !celltype %in% c("Unknown", "Total Cells in Sample")
  )
cell_type <- unique(scRNAseq_D28$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D28 <- list()
model_results_D28 <- list()


for (cell_type in unique(scRNAseq_D28$celltype)) {
  cell_type_data <- scRNAseq_D28 %>%
    filter(celltype == cell_type)

  D28_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D28[[cell_type]] <- D28_fit
  emm <- emmeans(D28_fit, ~ DPI)
  model_results_D28[[cell_type]] <- contrast(
    emm,
    method = list("D8-D2" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D2_D8.txt"
)

capture.output(
  lapply(model_results_D28, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.txt"
)

saveRDS(model_fit_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D2_D8.rds")
saveRDS(model_results_D28, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")

scRNAseq_D08 <- scRNAseq %>%
  filter(
    DPI %in% c("0 DPI", "8 DPI"),
   !celltype %in% c("Unknown", "Total Cells in Sample") 
  )
cell_type <- unique(scRNAseq_D08$celltype)
#Loop through each cell type at D0 and D2 and fit a linear mixed-effects model using lmer from the lme4 package
model_fit_D08 <- list()
model_results_D08 <- list()

for (cell_type in unique(scRNAseq_D08$celltype)) {
  cell_type_data <- scRNAseq_D08 %>%
    filter(celltype == cell_type)

  D08_fit <- lmer(
    Portion_of_PBMCs ~ DPI + (1 | Animal),
    data = cell_type_data
  )

  model_fit_D08[[cell_type]] <- D08_fit
  emm <- emmeans(D08_fit, ~ DPI)
  model_results_D08[[cell_type]] <- contrast(
    emm,
    method = list("D8-D0" = c(-1, 1)),
    adjust = "holm"
  )
}

capture.output(
  lapply(model_fit_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D8.txt"
)

capture.output(
  lapply(model_results_D08, summary),
  file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.txt"
)

saveRDS(model_fit_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_summaries_D0_D8.rds")
saveRDS(model_results_D08, file = "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
# Then make the plots for figures in the paper
#Will need to make celltype match in scRNAseq and Flow dataframes.
Flow <- read.csv("/scRNAseq/Sal_5pigs_2026/scSal_Flow/FlowProportons_CellConcentration_2025_02_13_scSalPaper_V2_Rready.csv", header = T)
#Make column "DPI" from "Sample"
Flow <- Flow %>%
  mutate(DPI = gsub(".*_(D\\d+)_.*", "\\1", Sample))
table(Flow$DPI)
# Set DPI as a factor
Flow$DPI <- as.factor(Flow$DPI)
#Make column "Animal" from "Sample", the 1st 3 digits of Sample are the animal number
Flow <- Flow %>%
  mutate(Animal = gsub("^(\\d{3}).*", "\\1", Sample))
#Make column "celltype" from "cell.type"
table(Flow$Animal)
# Set Animal as a factor
Flow$Animal <- as.factor(Flow$Animal)
head(Flow)
#Replace Flow$celltype==ASCs(Dump) with Flow$celltype==ASCs
Flow$celltype <- gsub("ASCs\\(Dump\\)", "ASCs", Flow$celltype)
table(Flow$celltype)
#For  Flow$celltype=="", replace with "CD3e- CD172a+ CD8a+"
Flow$celltype <- gsub("^$", "CD3e- CD172a+ CD8a+", Flow$celltype)
table(Flow$celltype)
#Change Flow$Percent_of_PBMCs name to Flow$Portion_of_PBMCs
Flow <- Flow %>%
  rename(Portion_of_PBMCs = Percent_of_PBMCs)

Flow$celltype <- gsub("CD2\\- GD T-cells", "CD2\\- gd T cells", Flow$celltype)
Flow$celltype <- gsub("CD4\\+ AB T-cells", "CD4\\+ ab T cells", Flow$celltype)
Flow$celltype <- gsub("monocytes CD4\\-", "Monocytes", Flow$celltype)
Flow$celltype <- gsub("pDCs & cDCs CD4\\+", "pDCs & cDCs", Flow$celltype)
Flow$celltype <- gsub("B-cells", "B cells", Flow$celltype)
table(Flow$celltype)

Flow_ready_merge <- Flow %>%
  filter( !celltype %in% c("CD8+ AB T-cells", "CD8+ CD4+ ab T cells","CD8- CD4- AB T-cells","CD2+ GD T-cells",  "CD3e- CD172a+ CD8a+","Neutrophils", "PBMCs") )
#add values for scRNAseq$celltype == pDCs & Flow$celltype == pDCs, then place into scRNAseq$celltype=="pDCs & cDCs" 
scRNAseq_nCells <- read.table("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Animal_Number_Prop_long.txt", sep = "\t", header = T)
scRNAseq_nCells$Sample_Celltype <- paste(scRNAseq_nCells$Sample, scRNAseq_nCells$Celltype, sep = "_")
head(scRNAseq_nCells)
scRNAseq_prop <- read.table("/scRNAseq/Sal_5pigs_2026/clustering_annotation/Block_cycling_genes/5pigs_salmonella_STACAS_via_Animal_ScranNorm_Regressed_PCA_d10_clustered_celltype_Sample100_Prop_long.txt", sep = "\t", header = T)
scRNAseq_prop$Sample_Celltype <- paste(scRNAseq_prop$Sample, scRNAseq_prop$Celltype, sep = "_")
head(scRNAseq_prop)
#merge the two dataframes by Sample_Celltype
scRNAseq <- merge(scRNAseq_nCells, scRNAseq_prop, by = "Sample_Celltype")
head(scRNAseq)
scRNAseq <- scRNAseq %>%
  rename(celltype = Celltype.x)
scRNAseq <- scRNAseq %>%
  rename(Sample = Sample.x)
scRNAseq$Celltype.y <- NULL
scRNAseq$Sample.y <- NULL
#Change scRNAseq$Proportion.of.Sample name to scRNAseq$Portion_of_PBMCs
scRNAseq <- scRNAseq %>%
  rename(Portion_of_PBMCs = Proportion.of.Sample)

scRNAseq$Technology <- "scRNAseq"
#scRNAseq$Sample values are formated as "###_D#".Pull the "D#" part to create the DPI column.
scRNAseq$DPI <- sub(".*_(D\\d+)$", "\\1", scRNAseq$Sample)
table(scRNAseq$DPI)
#scRNAseq$Sample values are formated as "###_D#".Pull the "###" part to create the Animal column.
scRNAseq$Animal <- gsub("^(\\d+)_D.*", "\\1", scRNAseq$Sample)
table(scRNAseq$Animal)
head(scRNAseq)
table(scRNAseq$celltype)
scRNAseq_pDC_cDC <- scRNAseq %>%
  filter(celltype %in% c("pDCs", "cDCs")) %>%
  group_by(Sample, DPI, Animal) %>%
  summarise(
    celltype = "pDCs & cDCs",
    nCells = sum(nCells, na.rm = TRUE),
    Portion_of_PBMCs = sum(Portion_of_PBMCs, na.rm = TRUE),
    .groups = "drop"
  )
scRNAseq_ready_merge <- scRNAseq %>%
  filter( !celltype %in% c("Unknown", "Total Cells in Sample") )%>%
  bind_rows(scRNAseq_pDC_cDC)
scRNAseq_ready_merge$Technology <- "scRNAseq"

table(Flow_ready_merge$celltype)
table(scRNAseq_ready_merge$celltype)
colnames(Flow_ready_merge)
colnames(scRNAseq_ready_merge)
table(Flow_ready_merge$Sample)
#Remove "_S" from the Flow_ready_merge$Sample values to match the scRNAseq_ready_merge$Sample values.
Flow_ready_merge$Sample <- gsub("_S", "", Flow_ready_merge$Sample)
table(Flow_ready_merge$Sample)
Flow_ready_merge$Sample_Celltype <- paste(Flow_ready_merge$Sample, Flow_ready_merge$celltype, sep = "_")
table(Flow_ready_merge$Sample_Celltype)
table(scRNAseq_ready_merge$Sample_Celltype)
#Trim whitespace from Flow_ready_merge$Sample_Celltype and scRNAseq_ready_merge$Sample_Celltype
Flow_ready_merge$Sample_Celltype <- trimws(as.character(Flow_ready_merge$Sample_Celltype))
scRNAseq_ready_merge$Sample_Celltype <- trimws(as.character(scRNAseq_ready_merge$Sample_Celltype))
colnames(Flow_ready_merge)
colnames(scRNAseq_ready_merge)
table(Flow_ready_merge$Technology)
# 105 rows
table(scRNAseq_ready_merge$Technology)
# 105 rows
#Add rows from Flow_ready_merge & scRNAseq_ready_merge to a new dataframe called Flow_scRNAseq
Flow_scRNAseq <- bind_rows(Flow_ready_merge, scRNAseq_ready_merge)
table(Flow_scRNAseq$Technology)
#flow cytometry       scRNAseq
#           105             105
# All transfered, now check to see if celltype names 
table(Flow_scRNAseq$celltype)
#change DPI format
Flow_scRNAseq$DPI <- gsub("D0", "0 DPI", Flow_scRNAseq$DPI)
Flow_scRNAseq$DPI <- gsub("D2", "2 DPI", Flow_scRNAseq$DPI)
Flow_scRNAseq$DPI <- gsub("D8", "8 DPI", Flow_scRNAseq$DPI)
table(Flow_scRNAseq$DPI)
#looks good save it
write.table(Flow_scRNAseq, "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_scRNAseq_merged_Portion_and_number_of_PBMCs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
Flow_scRNAseq <- read.table("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_scRNAseq_merged_Portion_and_number_of_PBMCs.txt", sep = "\t", header = TRUE)
Flow_scRNAseq$Percent_of_PBMCs <- Flow_scRNAseq$Portion_of_PBMCs * 100
#Add "%" to the end of the Flow_scRNAseq$Percent_of_PBMCs values
Flow_scRNAseq$Percent_of_PBMCs <- paste0(round(Flow_scRNAseq$Percent_of_PBMCs, 2), "%")
# make grouped bar plot of Portion_of_PBMCs for each celltype, grouped by DPI and colored by Technology
#Put all in 1 pdf to check that everything merged properly and that the celltype names match between Flow and scRNAseq
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_scRNAseq_merged_Portion_and_number_of_PBMCs.pdf", width = 12, height = 8)
ggplot(Flow_scRNAseq, aes(x = DPI, y = Portion_of_PBMCs, fill = Technology)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ celltype, scales = "free") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from  Flow cytometry and scRNAseq",
    x = "DPI",
    y = "Portion of PBMCs"
  ) +
  scale_fill_manual(values = c("flow cytometry" = "blue", "scRNAseq" = "red"))
dev.off()
#Everything looks good, now make grouped bar plot of Portion_of_PBMCs for Monocytes with error bars for each DPI and colored by Technology
Flow_scRNAseq_Monocytes <- Flow_scRNAseq %>%
  filter(celltype == "Monocytes")

pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_scRNAseq_merged_Portion_of_PBMC_Monocytes.pdf", width = 8, height = 6)
ggplot(Flow_scRNAseq_Monocytes, aes(x = DPI, y = Portion_of_PBMCs, fill = Technology)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "Portion of PBMCs; Cell-type distribution across 15 samples from  Flow cytometry and scRNAseq",
    x = "DPI",
    y = "Portion of PBMCs"
  ) 
dev.off()

Flow_scRNAseq_Monocytes_summary <- Flow_scRNAseq_Monocytes %>%
  group_by(DPI, Technology) %>%
  summarise(
    mean_value = mean(Portion_of_PBMCs, na.rm = TRUE),
    n = sum(!is.na(Portion_of_PBMCs)),
    se = sd(Portion_of_PBMCs, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
pdf("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_scRNAseq_merged_Portion_of_PBMC_Monocytes.pdf", width = 8, height = 6)
ggplot(Flow_scRNAseq_Monocytes_summary,aes(x = DPI, y = mean_value, fill = Technology)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  scale_fill_colorblind(black = FALSE) +
  theme_bw() +
  labs(
    title = "Portion of Monocytes in PBMCs",
    x = "DPI",
    y = "Portion of PBMCs",
    fill = "Technology"
  ) 
dev.off()

Flow_Only <- Flow_scRNAseq %>%
  filter(Technology == "flow cytometry")
cell_type <- unique(Flow_Only$celltype)
#Loop through each cell type 
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_boxplots.pdf",
  width = 6,
  height = 4
)
for (cell_type in unique(Flow_Only$celltype)) {
  cell_type_data <- Flow_Only %>%
    filter(celltype == cell_type)

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = Portion_of_PBMCs)
  ) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, outlier.shape = NA) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"), 
    ) +
    labs(
      title = paste(
        "Portion of PBMCs;",
        cell_type,
        "distribution across 15 samples from Flow cytometry"
      ),
      x = "DPI",
      y = "Portion of PBMCs"
    )

  print(plot)
}

dev.off()

pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_cells_per_ul_PBMCs_boxplots.pdf",
  width = 6,
  height = 4
)
for (cell_type in unique(Flow_Only$celltype)) {
  cell_type_data <- Flow_Only %>%
    filter(celltype == cell_type)

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = celltype_per_ul)
  ) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, outlier.shape = NA) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"), 
    ) +
    labs(
      title = paste(
        "Cells per ul;",
        cell_type,
        "distribution across 15 samples from Flow cytometry"
      ),
      x = "DPI",
      y = "Cells per ul"
    )

  print(plot)
}

dev.off()

scRNAseq_Only <- Flow_scRNAseq %>%
  filter(Technology == "scRNAseq")
cell_type <- unique(scRNAseq_Only$celltype)
#Loop through each cell type 
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_boxplots.pdf",
  width = 6,
  height = 4
)
for (cell_type in unique(scRNAseq_Only$celltype)) {
  cell_type_data <- scRNAseq_Only %>%
    filter(celltype == cell_type)

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = Portion_of_PBMCs)
  ) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, outlier.shape = NA) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"), 
    ) +
    labs(
      title = paste(
        "Portion of PBMCs;",
        cell_type,
        "distribution across 15 samples from scRNAseq"
      ),
      x = "DPI",
      y = "Portion of PBMCs"
    )

  print(plot)
}

dev.off()

Flow_nCells_D02_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")
p_values_Flow_nCells_D02 <- do.call(
  rbind,
  lapply(names(Flow_nCells_D02_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_nCells_D02_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

Flow_nCells_D28_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")
p_values_Flow_nCells_D28 <- do.call(
  rbind,
  lapply(names(Flow_nCells_D28_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_nCells_D28_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

Flow_nCells_D08_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
p_values_Flow_nCells_D08 <- do.call(
  rbind,
  lapply(names(Flow_nCells_D08_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_nCells_D08_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)
p_values_Flow_nCells_all <- rbind(p_values_Flow_nCells_D02, p_values_Flow_nCells_D28, p_values_Flow_nCells_D08)
write.table(p_values_Flow_nCells_all, "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_All_pvalues.txt", sep = "\t", row.names = FALSE, quote = FALSE)


Flow_prop_D02_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")
p_values_Flow_prop_D02 <- do.call(
  rbind,
  lapply(names(Flow_prop_D02_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_prop_D02_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

Flow_prop_D28_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")
p_values_Flow_prop_D28 <- do.call(
  rbind,
  lapply(names(Flow_prop_D28_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_prop_D28_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

Flow_prop_D08_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
p_values_Flow_prop_D08 <- do.call(
  rbind,
  lapply(names(Flow_prop_D08_emmeans_summary), function(cell_type) {
    result <- as.data.frame(Flow_prop_D08_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)
p_values_Flow_prop_all <- rbind(p_values_Flow_prop_D02, p_values_Flow_prop_D28, p_values_Flow_prop_D08)
write.table(p_values_Flow_prop_all, "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_All_pvalues.txt", sep = "\t", row.names = FALSE, quote = FALSE)


scRNAseq_prop_D02_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D2.rds")
p_values_scRNAseq_prop_D02 <- do.call(
  rbind,
  lapply(names(scRNAseq_prop_D02_emmeans_summary), function(cell_type) {
    result <- as.data.frame(scRNAseq_prop_D02_emmeans_summary[[cell_type]])

    data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

scRNAseq_prop_D28_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D2_D8.rds")
p_values_scRNAseq_prop_D28 <- do.call(
  rbind,
  lapply(names(scRNAseq_prop_D28_emmeans_summary), function(cell_type) {
    result <- as.data.frame(scRNAseq_prop_D28_emmeans_summary[[cell_type]])

  data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

scRNAseq_prop_D08_emmeans_summary <- readRDS("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_D0_D8.rds")
p_values_scRNAseq_prop_D08 <- do.call(
  rbind,
  lapply(names(scRNAseq_prop_D08_emmeans_summary), function(cell_type) {
    result <- as.data.frame(scRNAseq_prop_D08_emmeans_summary[[cell_type]])

  data.frame(
      CellType = cell_type,
      Contrast = result$contrast,
      estimate = result$estimate,
      SE = result$SE,
      df = result$df,
      t.ratio = result$t.ratio,
      p.value = result$p.value,
      row.names = NULL
    )
  })
)

p_values_scRNAseq_prop_all <- rbind(p_values_scRNAseq_prop_D02, p_values_scRNAseq_prop_D28, p_values_scRNAseq_prop_D08)
write.table(p_values_scRNAseq_prop_all, "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_All_pvalues.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# remove duplicate rows from "Flow_Only"
Flow_Only <- Flow_Only %>% distinct()
# remove duplicate rows from "p_values_Flow_prop_all"
p_values_Flow_prop_all <- p_values_Flow_prop_all %>% distinct()
pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_Portion_of_PBMCs_boxplots_with_significance.pdf",
  width = 6.2,
  height = 4
)
for (cell_type in unique(Flow_Only$celltype)) {
  cell_type_data <- Flow_Only %>%
    filter(celltype == cell_type)

  significance_data <- p_values_Flow_prop_all %>%
    filter(CellType == cell_type) %>%
    tidyr::separate(
      Contrast,
      into = c("group1", "group2"),
      sep = "\\s*-\\s*"
    ) %>%
    mutate(
      group1 = recode(
        group1,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      group2 = recode(
        group2,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      p = p.value,
      p.label = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*"
      ),
      y.position = max(
        cell_type_data$Portion_of_PBMCs,
        na.rm = TRUE
      ) * seq(1.10, 1.30, length.out = n())
    )

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = Portion_of_PBMCs)
  ) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, outlier.shape = NA) +
    ggpubr::stat_pvalue_manual(
      significance_data,
      label = "p.label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 7    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12 ,hjust = 0.5),
      plot.subtitle  = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    ) +
    labs(
      title = paste(
        "Portion of",
        cell_type,
        "in PBMCs from flow cytometry on whole blood"
      ),
      subtitle = "Significance: * p < 0.05   ** p < 0.01   *** p < 0.001",
      x = "DPI",
      y = "Portion of PBMCs \nfrom flow cytometry on whole blood"
    )

  print(plot)
}
dev.off()


p_values_Flow_nCells_all <- read.table("/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_Animal_Random_effect_model_emmeans_summaries_All_pvalues.txt", sep = "\t", header = TRUE)
# remove duplicate rows from "Flow_Only"
Flow_Only <- Flow_Only %>% distinct()
# remove duplicate rows from "p_values_Flow_nCells_all"
p_values_Flow_nCells_all <- p_values_Flow_nCells_all %>% distinct()

pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/Flow_celltype_per_ul_boxplots_with_significance.pdf",
  width = 6.2,
  height = 4
)
for (cell_type in unique(Flow_Only$celltype)) {
  cell_type_data <- Flow_Only %>%
    filter(celltype == cell_type)

  significance_data <- p_values_Flow_nCells_all %>%
    filter(CellType == cell_type) %>%
    tidyr::separate(
      Contrast,
      into = c("group1", "group2"),
      sep = "\\s*-\\s*"
    ) %>%
    mutate(
      group1 = recode(
        group1,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      group2 = recode(
        group2,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      p = p.value,
      p.label = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*"
      ),
      y.position = max(
        cell_type_data$celltype_per_ul,
        na.rm = TRUE
      ) * seq(1.10, 1.30, length.out = n())
    )

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = celltype_per_ul)
  ) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, outlier.shape = NA) +
    ggpubr::stat_pvalue_manual(
      significance_data,
      label = "p.label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 7    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle  = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    ) +
    labs(
      title = paste(
        "Cells per ul of",
        cell_type,
        "from flow cytometry on whole blood"
      ),
      subtitle = "Significance: * p < 0.05   ** p < 0.01   *** p < 0.001",
      x = "DPI",
      y = "Cells per ul \nfrom flow cytometry on whole blood"
    )

  print(plot)
}
dev.off()

p_values_scRNAseq_prop_all <- read.table("/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_Animal_Random_effect_model_emmeans_summaries_All_pvalues.txt", sep = "\t", header = TRUE)
# remove duplicate rows from "scRNAseq_Only"
scRNAseq_Only <- scRNAseq_Only %>% distinct()
# remove duplicate rows from "p_values_scRNAseq_prop_all"
p_values_scRNAseq_prop_all <- p_values_scRNAseq_prop_all %>% distinct()

pdf(
  "/scRNAseq/Sal_5pigs_2026/scSal_Flow/scRNAseq_Portion_of_PBMCs_boxplots_with_significance.pdf",
  width = 6.2,
  height = 4
)
for (cell_type in unique(scRNAseq_Only$celltype)) {
  cell_type_data <- scRNAseq_Only %>%
    filter(celltype == cell_type)

  significance_data <- p_values_scRNAseq_prop_all %>%
    filter(CellType == cell_type) %>%
    tidyr::separate(
      Contrast,
      into = c("group1", "group2"),
      sep = "\\s*-\\s*"
    ) %>%
    mutate(
      group1 = recode(
        group1,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      group2 = recode(
        group2,
        "D0" = "0 DPI",
        "D2" = "2 DPI",
        "D8" = "8 DPI"
      ),
      p = p.value,
      p.label = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*"
      ),
      y.position = max(
        cell_type_data$Portion_of_PBMCs,
        na.rm = TRUE
      ) * seq(1.10, 1.30, length.out = n())
    )

  plot <- ggplot(
    cell_type_data,
    aes(x = DPI, y = Portion_of_PBMCs)
  ) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, outlier.shape = NA) +
    ggpubr::stat_pvalue_manual(
      significance_data,
      label = "p.label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 7
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle  = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    ) +
    labs(
      title = paste(
        "Portion of",
        cell_type,
        "in PBMCs from scRNAseq"
      ),
      subtitle = "Significance: * p < 0.05   ** p < 0.01   *** p < 0.001",
      x = "DPI",
      y = "Portion of PBMCs"
    )

  print(plot)
}
dev.off()

sessionInfo()
# R version 4.5.3 (2026-03-11)
# Platform: x86_64-conda-linux-gnu
# Running under: AlmaLinux 9.6 (Sage Margay)

# Matrix products: default
# BLAS/LAPACK: /micromamba/envs/Emmeans/lib/libopenblasp-r0.3.34.so;  LAPACK version 3.12.0

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
#  [1] ggpubr_1.0.0      rstatix_1.1.0     gridExtra_2.3.1   hrbrthemes_0.9.3 
#  [5] viridis_0.6.5     viridisLite_0.4.3 ggthemes_6.0.0    ggplot2_4.0.3    
#  [9] dplyr_1.2.1       tidyr_1.3.2       emmeans_2.0.4     lme4_2.0-6       
# [13] Matrix_1.7-5     

# loaded via a namespace (and not attached):
#  [1] utf8_1.2.6         generics_0.1.4     stringi_1.8.9      lattice_0.23-1    
#  [5] magrittr_2.0.5     grid_4.5.3         estimability_2.0.0 RColorBrewer_1.1-3
#  [9] mvtnorm_1.4-2      backports_1.5.1    Formula_1.2-6      purrr_1.2.2       
# [13] scales_1.4.0       abind_1.4-8        reformulas_0.4.4   Rdpack_2.6.6      
# [17] cli_3.6.6          rlang_1.3.0        rbibutils_2.4.1    splines_4.5.3     
# [21] withr_3.0.3        parallel_4.5.3     pbkrtest_0.5.5     tools_4.5.3       
# [25] ggsignif_0.6.4     nloptr_2.2.1       minqa_1.2.8        boot_1.3-32       
# [29] broom_1.0.13       vctrs_0.7.3        R6_2.6.1           lifecycle_1.0.5   
# [33] stringr_1.6.0      car_3.1-5          MASS_7.3-66        pkgconfig_2.0.3   
# [37] pillar_1.11.1      gtable_0.3.6       glue_1.8.1         Rcpp_1.1.2        
# [41] tibble_3.3.1       tidyselect_1.2.1   farver_2.1.2       nlme_3.1-170      
# [45] labeling_0.4.3     carData_3.0-6      compiler_4.5.3     S7_0.2.2          