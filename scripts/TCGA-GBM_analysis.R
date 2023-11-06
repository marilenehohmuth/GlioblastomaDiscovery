# TCGA-GBM analysis #########################################################
# @ This script generates the plots present in the first main figure of the # 
# @ manuscript.                                                             #
#############################################################################


##########################
#### Loading packages ####
##########################

library(TCGAbiolinks)             # 2.23.5 -> 2.24.3
library(SummarizedExperiment)     # 1.20.0 -> 1.26.1
library(edgeR)                    # 3.32.1 -> 3.38.4
library(DESeq2)                   # ?      -> 1.36.0
library(ggplot2)                  # 3.3.5 
library(dplyr)                    # 1.0.7 
library(Matrix)                   # 1.3-4 
library(ggpubr)                   # 0.4.0
library(biomaRt)                  # 2.46.3 
library(GSVA)                     # 1.38.2
library(Hmisc)                    #
library(RColorBrewer)             #


# PART 1 - ANALYSIS OF GBM DATA FROM TCGA ######################################

################################
#### Step 1: Data retrieval ####
################################

# Retrieving data from the TCGA-GBM project.
query_GBM <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(
  query_GBM,
  method = "api",
  files.per.chunk = 10
)
data_GBM <- GDCprepare(query = query_GBM)

# Raw count data.
GBM_count_data <- as.data.frame(assay(data_GBM))

# Metadata.
GBM_metadata <- as.data.frame(colData(data_GBM))
  
#####################################################################
#### Step 2: Exploring PRNP expression in different tissue types ####
#####################################################################

# Normalizing raw count data of all samples, tumor and non-tumor (CPM normalization).
tissue_cpm <- cpm(GBM_count_data, log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
tissue_cpm1 <- tissue_cpm + 1

# Log-transforming CPM-normalized count data of all samples, tumor and non-tumor.
tissue_log <- log10(tissue_cpm1)

# Ensuring that the row names of the data.frame with log10(CPM+1) counts follow the same order
# as the row names of the data.frame with metadata.
tissue_log.t <- as.data.frame(t(tissue_log))
tissue_log.t <- tissue_log.t[rownames(GBM_metadata),]

# Mixing data.frames.
tissue_log_metadata <- cbind(GBM_metadata, tissue_log.t)

saveRDS(tissue_log_metadata, file = "results/TCGA/Tissue_logData+Metadata.RDS")

##########################################
#### Step 3: Selecting primary tumors ####
##########################################

# Filtering metadata.
GBM_primary_tumor_metadata <- GBM_metadata[GBM_metadata$definition == "Primary solid Tumor",]

# Removing samples of the neural subtype. 
GBM_primary_tumor_metadata <- GBM_primary_tumor_metadata %>% dplyr::mutate(
  paper_Transcriptome.Subtype_clean = ifelse(paper_Transcriptome.Subtype == 'NE', NA, as.character(paper_Transcriptome.Subtype))
)

# Filtering raw count data according to the metadata above.
GBM_primary_tumor_count_data <- GBM_count_data %>% dplyr::select(
  intersect(rownames(GBM_primary_tumor_metadata), colnames(GBM_count_data))
)

######################################################
#### Step 4: Normalizing primary tumor count data ####
######################################################

# Normalizing raw count data of primary tumors (CPM normalization).
pGBM_cpm <- cpm(GBM_primary_tumor_count_data, log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pGBM_cpm1 <- pGBM_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pGBM_log <- as.data.frame(log10(pGBM_cpm1))

#####################################################################################
#### Step 5: Exploring PRNP levels across different conditions in primary tumors ####
#####################################################################################

# Ensuring that the row names of the data.frame with log10(CPM+1) counts of primary tumors follow 
# the same order as the row names of the data.frame with metadata.
pGBM_log.t <- as.data.frame(t(pGBM_log))
pGBM_log.t <- pGBM_log.t[rownames(GBM_primary_tumor_metadata),]

# Mixing data.frames.
pGBM_log_metadata <- cbind(GBM_primary_tumor_metadata, pGBM_log.t)

saveRDS(pGBM_log_metadata, file = "results/TCGA/PrimaryGBMs_logData+Metadata.RDS")

#### @ FIGURE 1B, LEFT PANEL @ #### 
# Plotting PRNP levels across IDH-WT and IDH-mutant primary GBMs.
pGBM.idh.prnp <- ggplot(
  pGBM_log_metadata,
  aes(x = paper_IDH.status, y = ENSG00000171867.17, fill = paper_IDH.status)
) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  xlab(NA) +
  ggtitle("Primary GBMs") +
  ylab(expression(paste(italic("PRNP"), " expression"))) +
  scale_x_discrete(labels = c("IDH-mutant", "IDHwt", "Unclassified")) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#E5E0E4")) +
  stat_compare_means( # Kruskal-Wallis test
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "left",
    size = 4
  ) 

pdf("results/TCGA/PRNP_expression_GBM_Mutant_WT.pdf", width = 4, height = 4.5)
pGBM.idh.prnp
dev.off()

#### @ FIGURE 1B, RIGHT PANEL @ #### 
# Plotting PRNP levels across primary GBM subtypes.
pGBM.sub.prnp <- ggplot(
  pGBM_log_metadata,
  aes(x = paper_Transcriptome.Subtype_clean, y = ENSG00000171867.17, fill = paper_Transcriptome.Subtype_clean)
) + 
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  xlab(NA) +
  ylab(NA) + 
  stat_compare_means(
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "left",
    size = 4
  ) +
  scale_x_discrete(labels = c("Classical", "Mesenchymal", "Proneural", "Unclassified")) +
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#C77CFF", "#E5E0E4"))

pdf("results/TCGA/PRNP_expression_across_GBM_subtypes.pdf", width = 4, height = 4.5)
pGBM.sub.prnp
dev.off()

#### @ FIGURE 1B, LEFT+RIGHT PANEL @ #### 
pdf("results/TCGA/arranged_plots_pGBM.pdf", width = 8, height = 4)
cowplot::plot_grid(
  plotlist = list(pGBM.idh.prnp, pGBM.sub.prnp),
  align = "hv",
  ncol = 2,
  nrow = 1,
  axis = "btlr",
  rel_widths = c(1.8,2)
)
dev.off()

###############################################################
#### Step 6: Establishing groups with distinct PRNP levels ####
###############################################################

# Select row that contains log-normalized PRNP counts.
PRNP_counts_all <- pGBM_log["ENSG00000171867.17",]

# Transpose that row into a column and convert it to a dataframe.
PRNP_counts_all <- as.data.frame(t(PRNP_counts_all))
colnames(PRNP_counts_all) <- "PRNP"

# Calculate the quartiles of PRNP expression based on log-normalized counts.
PRNP_quartiles <- quantile(PRNP_counts_all$PRNP, c(0.25, 0.75))
PRNP_quartiles <- as.data.frame(PRNP_quartiles)
rownames(PRNP_quartiles) <- c("First", "Third")

# Classify samples according to their PRNP expression...
PRNP_counts_all <- PRNP_counts_all %>% dplyr::mutate(PRNP_status = case_when(
  PRNP <= PRNP_quartiles["First",] ~ "PRNP-Low", # Expression below the lower quartile = PRNP-Low
  PRNP > PRNP_quartiles["Third",] ~ "PRNP-High" # Expression above the upper quartile = PRNP-High
))

# Filtering out samples that do not belong to the PRNP-High and PRNP-Low groups.
PRNP_counts <- PRNP_counts_all[PRNP_counts_all$PRNP_status %in% c("PRNP-Low", "PRNP-High"),]

# This will be the raw count matrix for DESeq2.
PRNP_quartiles_count_data <- GBM_primary_tumor_count_data %>% dplyr::select(rownames(PRNP_counts))

# This will be the metadata for DESeq2.
PRNP_quartiles_status <- PRNP_counts %>% dplyr::select(PRNP_status)

##################################################################
#### Step 7: Characterizing the PRNP-High and PRNP-Low groups ####
##################################################################

# Selecting log10(CPM+1)-normalized counts of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_log_data <- as.data.frame(pGBM_log) %>% dplyr::select(rownames(PRNP_counts))

# Concatenating dataframes that contain the sample classification (PRNP-High/Low) and the log-normalized PRNP expression.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_status, t(PRNP_quartiles_log_data)) 

# Selecting metadata of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_metadata <- as.data.frame(t(GBM_primary_tumor_metadata)) %>% dplyr::select(rownames(PRNP_counts))

# Concatenating dataframes.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_log_metadata, t(PRNP_quartiles_metadata)) 
colnames(PRNP_quartiles_log_metadata)[colnames(PRNP_quartiles_log_metadata) == "ENSG00000171867.17"] <- "PRNP"

#### @ FIGURE 1C @ #### 
# Distribution of log-normalised PRNP counts in the PRNP-High/Low groups.
groups.prnp.exp.hist <- ggplot(
  PRNP_counts,
  aes(x = PRNP, fill = PRNP_status)
) +
  geom_histogram(color = "black", alpha=0.6, position = "identity", binwidth = 0.025) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  ylab("Frequency") +
  xlab("PRNP expression") +
  scale_fill_manual(name = "Group", values = c("orange", "steelblue1")) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )

pdf("results/TCGA/PRNP_distribution_PRNP-High_and_PRNP-Low_groups.pdf", width = 7, height = 3)
groups.prnp.exp.hist
dev.off()

#### @ FIGURE 1D, LEFT PANEL @ #### 
# Proportion of IDH-WT, IDH-mutant and unclassified samples in the PRNP-High/Low groups.
groups.prnp.idh <- ggplot(
  PRNP_quartiles_log_metadata,
  aes(x = reorder(PRNP_status, PRNP), fill = unlist(paper_IDH.status))
) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  xlab(NA) +
  ylab("Proportion") +
  theme(
    axis.text.x = element_text(angle = 0, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  scale_fill_manual(
    name = "IDH status",
    values = c("#00BFC4", "#F8766D", "snow2"),
    labels = c("IDH-mutant", "IDHwt", "Unclassified")
  ) +
  scale_x_discrete(labels = c(expression(italic("PRNP")^"low"), expression(italic("PRNP")^"high")))

pdf("results/TCGA/IDH_status_across_PRNP-High_and_PRNP-Low_groups.pdf", width = 4, height = 4.5)
groups.prnp.idh
dev.off()

#### @ FIGURE 1D, RIGHT PANEL @ #### 
# Proportion of classical, mesenchymal, proneural and unclassified samples in the PRNP-High/Low groups.
groups.prnp.sub <- ggplot(
  PRNP_quartiles_log_metadata,
  aes(x = reorder(PRNP_status, PRNP), fill = unlist(paper_Transcriptome.Subtype_clean))
) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  xlab(NA) +
  ylab("Proportion") + 
  theme(
    axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  scale_fill_manual(
    name = "Subtype", 
    values = c("#F8766D", "#7CAE00", "#C77CFF", "snow2"),
    labels = c("Classical","Mesenchymal", "Proneural", "Unclassified")
  ) +
  scale_x_discrete(labels = c(expression(italic("PRNP")^"low"), expression(italic("PRNP")^"high")))

pdf("results/TCGA/Subtypes_across_PRNP-High_and_PRNP-Low_groups.pdf", width = 4, height = 4.5)
groups.prnp.sub
dev.off()

##################################################################
#### Step 8: Identifying differentially expressed transcripts ####
##################################################################

# Creating DESeq2 object.
dds <- DESeqDataSetFromMatrix(
  countData = PRNP_quartiles_count_data, 
  colData = PRNP_quartiles_status,
  design = ~ PRNP_status
)

# Filtering out genes that are expressed in very low levels.
dds <- dds[rowSums(counts(dds)) > 1,]

# Estimating size factors & running DESeq2 analysis.
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Extracting results with contrast of interest.
res_unshrunken <- results(
  dds, 
  contrast = c('PRNP_status','PRNP-High', 'PRNP-Low'), 
  alpha = 0.05
)

# Doing log fold change shrinkage based on a normal distribution.
res_shrunken <- lfcShrink(
  dds = dds, 
  res = res_unshrunken, 
  contrast = c('PRNP_status','PRNP-High', 'PRNP-Low'), 
  type = "normal"
)

# Keeping only statistically significant results.
res_shrunken <- as.data.frame(res_shrunken)
deg.sig <- res_shrunken[res_shrunken$padj <= 0.05,]
write.csv(deg.sig, file = "results/TCGA/TCGA-GBM_DETs_padj005.csv")

# Adding gene classification to complete dataframe.
deg.all <- res_shrunken %>% mutate(classification = case_when(
  padj <= 0.05 & log2FoldChange > 0 ~ "Upregulated",
  padj <= 0.05 & log2FoldChange < 0 ~ "Downregulated",
  padj > 0.05 ~ "Non-significant",
  is.na(padj) ~ "Non-significant")
)

labs <- c(
  paste0("Upregulated (n=", nrow(deg.all[deg.all$classification == "Upregulated",]), ")"),
  paste0("Non-significant (n=", nrow(deg.all[deg.all$classification == "Non-significant",]), ")"),
  paste0("Downregulated (n=", nrow(deg.all[deg.all$classification == "Downregulated",]), ")")
)

#### @ FIGURE 1E @ ####
# Volcano plot with differentially expressed transcripts. 
volcano <- ggplot(
  deg.all,
  aes(x = log2FoldChange, y = -log10(padj), color = classification)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_manual(
    values = c("blue", "lightgray", "red"),
    name = "",
    labels = labs) +
  xlab(expression(log[2]("Fold change"))) +
  ylab(expression(-log[10]("Adjusted p-value"))) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.text.align = 0
  ) + 
  guides(color = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

pdf("results/TCGA/volcano.pdf", width = 6, height = 4)
volcano
dev.off()

#####################################################
#### Step 9: Doing over-representation analysis  ####
#####################################################

# Need to add over-representation analysis or GSEA with clusterProfiler inside R to remove external dependency on 
# g:Profiler (this is a website where you can upload a gene list and get enrichment analysis results).

#### @ FIGURE 1F @ ####
# Will be a bar plot showing intracellular traffic and vesicle-related terms.