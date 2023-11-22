# TCGA-GBM ANALYSIS #########################################################
# @ This script generates the plots present in the first main figure of the # 
# @ manuscript.                                                             #
#############################################################################


##########################
#### Loading packages ####
##########################

# R 4.3.2
library(TCGAbiolinks)             # 2.28.3
library(SummarizedExperiment)     # 1.30.2
library(edgeR)                    # 3.42.4
library(DESeq2)                   # 1.40.2
library(ggplot2)                  # 3.4.4
library(dplyr)                    # 1.1.3
library(ggpubr)                   # 0.6.0
library(biomaRt)                  # 2.56.1
library(GSVA)                     # 1.48.2
library(Hmisc)                    # 5.1_1
library(RColorBrewer)             # 1.1_3
library(clusterProfiler)          # 4.8.1
library(org.Hs.eg.db)             # 3.17.0
library(ggrepel)                  # 0.9.4

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

saveRDS(tissue_log_metadata, file = paste0(getwd(), "/results/TCGA-GBM/Tissue_logData+Metadata.RDS"))

##########################################
#### Step 3: Selecting primary tumors ####
##########################################

# Filtering metadata.
GBM_primary_tumor_metadata <- GBM_metadata[GBM_metadata$definition == "Primary solid Tumor",]

# Creating new subtype column where the neural subtype will be classified as NA (unknown subtype).
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
##################################################################################################
#                                                                                                #
# Note: at this point, we will keep samples whose IDH status is mutant/unknown and samples whose #
# transcriptional subtype is unknown (unknown now also includes the neural samples!) just to     #
# to see the expression patterns of PRNP in these conditions. However, these samples will be     # 
# excluded from the differential gene expression analysis that will be performed later on.       #                           
#                                                                                                #
##################################################################################################

# Ensuring that the row names of the data.frame with log10(CPM+1) counts of primary tumors follow 
# the same order as the row names of the data.frame with metadata.
pGBM_log.t <- as.data.frame(t(pGBM_log))
pGBM_log.t <- pGBM_log.t[rownames(GBM_primary_tumor_metadata),]

# Mixing data.frames.
pGBM_log_metadata <- cbind(GBM_primary_tumor_metadata, pGBM_log.t)

saveRDS(pGBM_log_metadata, file = paste0(getwd(), "/results/TCGA-GBM/PrimaryGBMs_logData+Metadata.RDS"))

#### @ FIGURE 1B, LEFT PANEL (MAIN) @ #### 
# Plotting PRNP levels across IDH-WT and IDH-mutant primary GBMs.
pGBM.idh.prnp <- ggplot(
  pGBM_log_metadata,
  aes(x = paper_IDH.status, y = ENSG00000171867.17)
) +
  geom_violin(aes(fill = paper_IDH.status), scale = "width") +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  ylab(expression("Normalized"~italic("PRNP")~"expression")) +
  scale_fill_manual(
    name = "IDH status",
    values = c("#00BFC4", "#F8766D", "#E5E0E4"),
    labels = c(
      paste0("IDH-mutant (n=", nrow(pGBM_log_metadata[pGBM_log_metadata$paper_IDH.status == "Mutant",]), ")"), 
      paste0("IDHwt (n=", nrow(pGBM_log_metadata[pGBM_log_metadata$paper_IDH.status == "WT",]), ")"),
      paste0("Unclassified (n=", nrow(pGBM_log_metadata[is.na(pGBM_log_metadata$paper_IDH.status),]), ")")
    )
  ) +
  stat_compare_means( # Kruskal-Wallis test
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "left",
    size = 4
  ) 

pdf(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_primaryTumours_PRNP_expression_perIDHstatus.pdf"), 
  width = 5, 
  height = 3
)
pGBM.idh.prnp
dev.off()

#### @ FIGURE 1B, RIGHT PANEL (MAIN) @ #### 
# Plotting PRNP levels across primary GBM subtypes.
pGBM.sub.prnp <- ggplot(
  pGBM_log_metadata,
  aes(x = paper_Transcriptome.Subtype_clean, y = ENSG00000171867.17)
) + 
  geom_violin(aes(fill = paper_Transcriptome.Subtype_clean), scale = "width") +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  xlab(NULL) +
  ylab(expression("Normalized"~italic("PRNP")~"expression")) + 
  stat_compare_means( # Kruskal-Wallis test
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "left",
    size = 4
  ) +
  scale_fill_manual(
    name = "Transcriptional subtype",
    values = c("#F8766D", "#7CAE00", "#C77CFF", "#E5E0E4"),
    labels = c(
      paste0("Classical (n=", nrow(pGBM_log_metadata[pGBM_log_metadata$paper_Transcriptome.Subtype_clean == "CL",]), ")"), 
      paste0("Mesenchymal (n=", nrow(pGBM_log_metadata[pGBM_log_metadata$paper_Transcriptome.Subtype_clean == "ME",]), ")"),
      paste0("Proneural (n=", nrow(pGBM_log_metadata[pGBM_log_metadata$paper_Transcriptome.Subtype_clean == "PN",]), ")"),
      paste0("Unclassified (n=", nrow(pGBM_log_metadata[is.na(pGBM_log_metadata$paper_Transcriptome.Subtype_clean),]), ")")
    )
  )

pdf(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_primaryTumours_PRNP_expression_perTranscriptionalSubtype.pdf"), 
  width = 5, 
  height = 3
)
pGBM.sub.prnp
dev.off()

#### @ FIGURE 1B, LEFT+RIGHT PANEL (MAIN) @ #### 
pdf(paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_primaryTumours_PRNP_expression_perIDHstatus_plus_perSubtype_arranged_plots.pdf"), width = 10, height = 3)
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
##################################################################################################
#                                                                                                #
# Now that we explored PRNP expression levels in all IDH and transcriptional subtype conditions, #
# we will exclude samples whose IDH status is mutant/unknown and samples whose transcriptional   #
# subtype is unknown. After that, we will create groups of samples that express high levels of   # 
# PRNP (PRNP-High) or low levels of PRNP (PRNP-Low) and explore those groups a little bit.       #                           
#                                                                                                #
##################################################################################################

# Remove samples whose IDH status is mutant/unknown.
pGBM_metadata_filt <- GBM_primary_tumor_metadata %>% dplyr::filter(paper_IDH.status == "WT")

# Remove samples whose transcriptional subtype is unknown.
pGBM_metadata_filt <- pGBM_metadata_filt %>% dplyr::filter(!is.na(paper_Transcriptome.Subtype_clean))

# Filter count data based on filtered metadata.
pGBM_count_data_filt <- GBM_primary_tumor_count_data %>% dplyr::select(rownames(pGBM_metadata_filt))

# Normalize filtered count data according to log10(CPM+1).
pGBM_filt_norm <- cpm(pGBM_count_data_filt, log = FALSE)
pGBM_filt_norm <- pGBM_filt_norm + 1
pGBM_filt_norm <- log10(pGBM_filt_norm)

# Select row that contains log-normalized PRNP counts.
PRNP_counts_all <- pGBM_filt_norm["ENSG00000171867.17",]

# Put the numeric vector with log-normalized PRNP counts into a dataframe.
PRNP_counts_all <- data.frame(PRNP = PRNP_counts_all)

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
PRNP_quartiles_count_data <- pGBM_count_data_filt %>% dplyr::select(rownames(PRNP_counts))

# This will be the metadata for DESeq2.
PRNP_quartiles_status <- PRNP_counts %>% dplyr::select(PRNP_status)

##################################################################
#### Step 7: Characterizing the PRNP-High and PRNP-Low groups ####
##################################################################

# Selecting log10(CPM+1)-normalized counts of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_log_data <- as.data.frame(pGBM_filt_norm) %>% dplyr::select(rownames(PRNP_counts))

# Concatenating dataframes that contain the sample classification (PRNP-High/Low) and the log-normalized PRNP expression.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_status, t(PRNP_quartiles_log_data)) 

# Selecting metadata of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_metadata <- as.data.frame(t(pGBM_metadata_filt)) %>% dplyr::select(rownames(PRNP_counts))

# Concatenating dataframes.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_log_metadata, t(PRNP_quartiles_metadata)) 
colnames(PRNP_quartiles_log_metadata)[colnames(PRNP_quartiles_log_metadata) == "ENSG00000171867.17"] <- "PRNP"

#### @ FIGURE 1C (MAIN) @ #### 
# Distribution of log-normalised PRNP counts in the PRNP-High/Low groups.
groups.prnp.exp.hist <- ggplot(
  PRNP_counts,
  aes(x = PRNP, fill = PRNP_status)
) +
  geom_histogram(color = "black", alpha=0.6, position = "identity", binwidth = 0.025) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  xlab(expression(paste(italic("PRNP"), " expression"))) +
  ylab("Frequency") +
  scale_fill_manual(
    name = "Group", 
    values = c("orange", "steelblue1"),
    labels = c(expression(italic("PRNP")^"high"), expression(italic("PRNP")^"low"))) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.text.align = 0
  )

pdf(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_PRNP_expression_histogram.pdf"), 
  width = 7, 
  height = 3
)
groups.prnp.exp.hist
dev.off()

#### @ FIGURE 1D, RIGHT PANEL (MAIN) @ #### 
# Proportion of classical, mesenchymal, proneural and unclassified samples in the PRNP-High/Low groups.
groups.prnp.sub <- ggplot(
  PRNP_quartiles_log_metadata,
  aes(x = reorder(PRNP_status, PRNP), fill = unlist(paper_Transcriptome.Subtype_clean))
) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  xlab(NULL) +
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

pdf(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_composition_subtypes.pdf"), 
  width = 4, 
  height = 4.5
)
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

deg.all <- as.data.frame(res_shrunken)
deg.all <- deg.all %>% dplyr::filter(!is.na(padj))

# Creating column with ENSEMBL IDs (without version identifier).
deg.all$ensembl <- gsub("\\..*", "", rownames(deg.all))

# Getting correspondence between ENSEMBL IDs and gene symbols.
correspondence <- select(
  org.Hs.eg.db, 
  keys = deg.all$ensembl, 
  columns = "SYMBOL", 
  keytype = "ENSEMBL"
)

# Adding column with gene symbols.
deg.all$gene_symbol <- correspondence$SYMBOL[match(deg.all$ensembl, correspondence$ENSEMBL)]
deg.all$gene_symbol[is.na(deg.all$gene_symbol)] <- "Unknown"

# Adding gene classification to complete dataframe.
deg.all <- deg.all %>% mutate(classification = case_when(
  padj <= 0.05 & log2FoldChange > 0 ~ "Upregulated",
  padj <= 0.05 & log2FoldChange < 0 ~ "Downregulated",
  padj > 0.05 ~ "Non-significant",
  is.na(padj) ~ "Non-significant")
)

# Saving only statistically significant results.
deg.sig <- deg.all %>% dplyr::filter(padj <= 0.05)
write.csv(deg.sig, file = paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_DETs_padj0.05.csv"))

labs <- c(
  paste0("Downregulated (n=", nrow(deg.all[deg.all$classification == "Downregulated",]), ")"),
  paste0("Non-significant (n=", nrow(deg.all[deg.all$classification == "Non-significant",]), ")"),
  paste0("Upregulated (n=", nrow(deg.all[deg.all$classification == "Upregulated",]), ")")
)

#### @ FIGURE 1E (MAIN) @ ####
# Volcano plot with differentially expressed transcripts. 
volcano <- ggplot(
  deg.all[deg.all$gene_symbol != "PRNP",],
  aes(x = log2FoldChange, y = -log10(padj), color = classification)) +
  geom_point(size = 1, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(
    values = c("royalblue1", "lightgray", "indianred2"),
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
    legend.text.align = 0,
    plot.title = element_text(size = 15, hjust = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggtitle(expression(italic("PRNP")^"high"~italic("versus")~italic("PRNP")^"low"))

pdf(paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_volcano_withoutPRNP.pdf"), width = 8, height = 4)
volcano
dev.off()

# Also create an alternative plot in which PRNP is not removed.
volcano <- ggplot(
  deg.all,
  aes(x = log2FoldChange, y = -log10(padj), color = classification)) +
  geom_point(size = 1, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(
    values = c("royalblue1", "lightgray", "indianred2"),
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
    legend.text.align = 0,
    plot.title = element_text(size = 15, hjust = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggtitle(expression(italic("PRNP")^"high"~italic("versus")~italic("PRNP")^"low")) +
  geom_text_repel(
    data = deg.all[deg.all$gene_symbol == "PRNP",],
    aes(label = gene_symbol),
    color = "black"
  )

pdf(paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_volcano_withPRNP.pdf"), width = 8, height = 4)
volcano
dev.off()

#######################################
#### Step 9: Functional profiling  ####
#######################################

#### GSEA ####

# Get ranked gene list.
gene_list <- deg.sig$log2FoldChange
names(gene_list) <- deg.sig$ensembl
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

# Perform Gene Set Enrichment Analysis (GSEA).
gse <- gseGO(
  geneList = gene_list, 
  ont = "ALL", 
  keyType = "ENSEMBL", 
  nPerm = 10000, 
  minGSSize = 3, 
  maxGSSize = 800, 
  pvalueCutoff = 0.05, 
  verbose = TRUE, 
  OrgDb = org.Hs.eg.db, 
  pAdjustMethod = "BH"
)
gse <- as.data.frame(gse@result)

# Save GSEA results to output file.
write.csv(
  gse, 
  file = paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_GSEA_padj0.05_BH_correctionMethod.csv")
)

terms <- c(
  "immunoglobulin complex",
  "immune response",
  "motile cilium",
  "ciliary plasm",
  "external side of plasma membrane",
  "cilium movement",
  "microtubule bundle formation",
  "plasma membrane bounded cell projection assembly",
  "cell projection assembly",
  "positive regulation of epidermal growth factor receptor signaling pathway",
  "inner dynein arm assembly",
  "positive regulation of epidermal growth factor-activated receptor activity",
  "extracellular matrix",
  "extracellular region",
  "synaptic vesicle membrane",
  "exocytic vesicle membrane",
  "cell periphery",
  "plasma membrane",
  "double-strand break repair via homologous recombination",
  "regulation of cell cycle phase transition",
  "double-strand break repair",
  "DNA damage response",
  "cell division",
  "cell cycle checkpoint signaling",
  "negative regulation of cell cycle process",
  "DNA repair",
  "neuron differentiation",
  "chromosome organization"
)

#### @ FIGURE 1F (MAIN) @ ####
pdf(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_GSEA_padj0.05_BH_correctionMethod_SelectedTerms.pdf"), 
  width = 10, 
  height = 8
)
ggplot(
  gse[gse$Description %in% terms,],
  aes(x = NES, y = reorder(Description, NES), fill = -log10(p.adjust), size = setSize)
) +
  geom_point(shape = 21, color = "black") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  xlab("Normalized enrichment scores (NES)") +
  ylab("Gene Ontology (GO) term") +
  scale_fill_gradientn(
    colors = colorRampPalette(brewer.pal(5, "Reds"))(100),
    name = expression(-log[10]~"(Adjusted p-value)")
  ) +
  scale_size_continuous(range = c(1,5), name = "# Genes") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle(expression("GSEA -"~italic("PRNP")^"high"~italic("versus")~italic("PRNP")^"low"))
dev.off()

#### ORA ####

# Perform Over-Representation Analysis (ORA) on upregulated transcripts.
ego_up <- enrichGO(
  gene = names(gene_list)[gene_list > 0],
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = "ENSEMBL",
  readable = TRUE
)

# Save ORA results to output file.
write.csv(
  as.data.frame(ego_up@result), 
  file = paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_ORA_upregulated_transcripts_padj0.05_BH_correctionMethod.csv")
)

# Perform Over-Representation Analysis (ORA) on downregulated transcripts.
ego_down <- enrichGO(
  gene = names(gene_list)[gene_list < 0],
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = "ENSEMBL",
  readable = TRUE
)

# Save ORA results to output file.
write.csv(
  as.data.frame(ego_down@result), 
  file = paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_ORA_downregulated_transcripts_padj0.05_BH_correctionMethod.csv")
)