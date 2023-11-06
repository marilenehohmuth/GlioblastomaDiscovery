# Analysis of single-cell datasets ###########################################
# @ This script generates the plots present in the second main figure of the # 
# @ manuscript + supplemental material.                                      #
##############################################################################


##########################
#### Loading packages ####
##########################

library(fusca)                    # 1.2.1
library(ggplot2)                  # 3.3.5
library(dplyr)                    # 1.0.7
library(Matrix)                   # 1.3-4
library(RColorBrewer)             # 1.1-2
library(Hmisc)                    # 4.5-0
library(corrplot)                 # 0.92

###################################
###################################
#---- Darmanis et al. dataset ----#
###################################
###################################

##########################################################
#### Step 1: Loading and subsetting data for analysis ####
##########################################################

# Load Darmanis et al count data.
data_darmanis <- read.csv(
    "data/darmanis/GBM_raw_gene_counts.csv",
    header = TRUE,
    row.names = 1,
    sep = " ",
    check.names = FALSE
)

# Load Darmanis et al metadata.
metadata_darmanis <- read.csv(
    "data/darmanis/GBM_metadata.csv",
    header = TRUE,
    row.names = 1,
    sep = " "
)

# Subset count data to keep only malignant cells.
malignant_data_darmanis <- data_darmanis[,colnames(data_darmanis) %in% 
    rownames(metadata_darmanis)[metadata_darmanis$Cluster_2d %in% c(1, 4, 11)]]

# Subset metadata to keep only malignant cells.
malignant_metadata_darmanis <- metadata_darmanis[metadata_darmanis$Cluster_2d %in% c(1, 4, 11),]

# Convert count data data.frame to sparse matrix.
malignant_data_sparse_darmanis <- as(as.matrix(malignant_data_darmanis), "sparseMatrix")

#################################
#### Step 2: Data processing ####
#################################

# Create CellRouter object for Darmanis et al dataset.
cellrouter_darmanis <- CreateCellRouter(
    malignant_data_sparse_darmanis,
    min.genes = 0,
    min.cells = 0,
    is.expr = 0
)
cellrouter_darmanis@assays$RNA@rawdata <- cellrouter_darmanis@assays$RNA@rawdata[rownames(cellrouter_darmanis@assays$RNA@ndata), colnames(cellrouter_darmanis@assays$RNA@ndata)]

# Ensure row names in metadata table have the same order as row names in the sampTab slot of the RNA assay & add
# sample identification to the CellRouter object.
malignant_metadata_darmanis <- malignant_metadata_darmanis[rownames(cellrouter_darmanis@assays$RNA@sampTab),]
cellrouter_darmanis@assays$RNA@sampTab$samples <- malignant_metadata_darmanis$Sample.name

# Normalize malignant count data.
cellrouter_darmanis <- Normalize(cellrouter_darmanis)

# Transpose normalized data & ensure row names have the same order as row names in the sampTab slot of the RNA assay.
transposed_norm_data <- as.data.frame(t(as.matrix(cellrouter_darmanis@assays$RNA@ndata)))
transposed_norm_data <- transposed_norm_data[rownames(cellrouter_darmanis@assays$RNA@sampTab),]

# Concatenate dataframes to have gene expression and metadata information in the same table.
expr <- cbind(t(as.matrix(cellrouter_darmanis@assays$RNA@ndata)), cellrouter_darmanis@assays$RNA@sampTab)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
boxplot.samples.darmanis <- ggplot(
    expr,
    aes(x = reorder(samples, PRNP), y = PRNP, fill = samples)
) +
    geom_boxplot(color = "antiquewhite4", outlier.colour = "antiquewhite4") +
    xlab("Sample") +
    ylab("Normalized PRNP counts") +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1)
    )

# Scale data according to all genes.
cellrouter_darmanis <- scaleData(cellrouter_darmanis, blocksize = nrow(cellrouter_darmanis@assays$RNA@ndata))

# Perform PCA.
cellrouter_darmanis <- computePCA(cellrouter_darmanis, num.pcs = 50, seed = 42)

# Compute UMAP.
cellrouter_darmanis <- computeUMAP(cellrouter_darmanis, num.pcs = 15, seed = 42)

# Save processed CellRouter object of the Darmanis et al dataset.
save(cellrouter_darmanis, file = "results/darmanis/Darmanis_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")

###############################################################
#### Step 3: Establishing groups with distinct PRNP levels ####
###############################################################

# Getting cells that express PRNP.
prnp_pos_darmanis <- colnames(cellrouter_darmanis@assays$RNA@ndata)[cellrouter_darmanis@assays$RNA@ndata["PRNP",] > 0]
# Getting cells that do not express PRNP.
prnp_neg_darmanis <- colnames(cellrouter_darmanis@assays$RNA@ndata)[cellrouter_darmanis@assays$RNA@ndata["PRNP",] == 0]

cellrouter_darmanis@assays$RNA@sampTab <- cellrouter_darmanis@assays$RNA@sampTab %>% mutate(PRNP_status = case_when(
    rownames(cellrouter_darmanis@assays$RNA@sampTab) %in% prnp_pos_darmanis ~ "PRNP_positive_cells",
    rownames(cellrouter_darmanis@assays$RNA@sampTab) %in% prnp_neg_darmanis ~ "PRNP_negative_cells"
))

# Ensure row names in cell embeddings of UMAP slot have the same order as the row names of the transposed, normalized count data.
cellrouter_darmanis@umap$cell.embeddings <- cellrouter_darmanis@umap$cell.embeddings[rownames(transposed_norm_data),]

# Concatenate dataframes to have cell embeddings and gene expression information in the same table.
dim.red.darmanis <- cbind(cellrouter_darmanis@umap$cell.embeddings, transposed_norm_data)
colnames(dim.red.darmanis)[1] <- "UMAP1"
colnames(dim.red.darmanis)[2] <- "UMAP2"

# Ensure row names have the same order as the row names in the sampTab slot of the RNA assay.
dim.red.darmanis <- dim.red.darmanis[rownames(cellrouter_darmanis@assays$RNA@sampTab),]

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.darmanis <- cbind(dim.red.darmanis, cellrouter_darmanis@assays$RNA@sampTab)

#### @ FIGURE 2B, LEFT PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
umap_darmanis_prnp <- ggplot(
    dim.red.darmanis %>% dplyr::arrange(PRNP),
    aes(x = UMAP1, y = UMAP2, color = PRNP)
) +
    geom_point(size = 1) +
    theme_classic() +
    theme(
        legend.position = "right",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)
    ) +
    scale_color_gradientn(
        colors = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(1000),
        limits = c(min(dim.red.darmanis$PRNP), max(dim.red.darmanis$PRNP))
    )

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample ID.
umap_darmanis_sample <- ggplot(
    dim.red.darmanis,
    aes(x = UMAP1, y = UMAP2, color = samples)
) +
    geom_point(size = 1) +
    theme_void() +
    theme(
        legend.position = "bottom",
        legend.key.size = unit(2, 'cm'),
        legend.box.just = "center"
    ) + 
    guides(
        color = guide_legend(override.aes = list(size=3),
        keywidth = 0,
        keyheight = 0,
        ncol = 4,
        title = "")
    ) 

#############################################################
#### Step 4: Finding signatures of PRNP+ and PRNP- cells ####
#############################################################

# Find PRNP+ and PRNP- marker gene signatures & save results to output file.
markers_darmanis <- findSignatures(
    object = cellrouter_darmanis,
    assay.type = "RNA",
    column = "PRNP_status",
    test.use = "wilcox"
)
write.csv(markers_darmanis, "results/darmanis/Darmanis_PRNP+_vs_PRNP-_signature_markers.csv")

# This is a volcano plot that, at the moment, is not shown in the paper.
ggplot(
    markers_darmanis,
    aes(x = fc, y = -log10(pval))) +
  geom_point(size = 1, fill = "royalblue", shape = 21) +
  theme_classic() + 
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
    ) +
  scale_y_break(c(50,185), ticklabels=c(190)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ylab(expression(-log[10](p-value))) +
  xlab(expression(log[2]("Fold change")))

######################################
#### Step 5: Correlation analysis ####
######################################

df_darmanis <- data.frame(matrix(0, ncol=4))
colnames(df_darmanis) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_darmanis@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_darmanis@assays$RNA@ndata[g,]), as.numeric(cellrouter_darmanis@assays$RNA@ndata['PRNP',]))
  df_darmanis[i,'gene1'] <- 'PRNP'
  df_darmanis[i,'gene2'] <- g
  df_darmanis[i,'correlation'] <- c$estimate # Store the correlation.
  df_darmanis[i,'pvalue'] <- c$p.value # Stores the p-value.
  i <- i + 1
}
write.csv(df_darmanis, file = "results/darmanis/Darmanis_Correlation_PRNPvsAllGenes.csv")

df_darmanis <- df_darmanis[!is.na(df_darmanis$correlation),]
df_darmanis <- df_darmanis[order(df_darmanis$correlation, decreasing = TRUE),]
df_darmanis <- df_darmanis[-c(1),]
rownames(df_darmanis) <- 1:length(rownames(df_darmanis))

df_darmanis <- df_darmanis %>% mutate(significant = case_when(
    pvalue <= 0.05 ~ "Significant",
    pvalue > 0.05 ~ "Non-significant"
))

#### @ FIGURE 2C, LEFT PANEL (MAIN) @ #### 
darmanis.cor <- ggplot(
    df_darmanis,
    aes(y = correlation, x = 1:length(rownames(df_darmanis)), color = factor(significant))
) +
    geom_point() +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ylab("Pearson correlation") +
    xlab("Number of genes") +
    scale_color_manual(values = c("gray", "deeppink"), name = "Significance")  +
    ggtitle("Darmanis et al") +
    theme(plot.title = element_text(hjust = 0.5))

#################################
#################################
#---- Neftel et al. dataset ----#
#################################
#################################

##########################################################
#### Step 1: Loading and subsetting data for analysis ####
##########################################################

# Load Neftel et al count data.
data_neftel <- read.table(
    "data/neftel/IDHwtGBM.processed.SS2.logTPM.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

# Load Neftel et al metadata.
metadata_neftel <- read.table(
    "data/neftel/IDHwt.GBM.Metadata.SS2.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

# Subset count data to keep only malignant cells.
malignant_data_neftel <- data_neftel[,colnames(data_neftel) %in% rownames(metadata_neftel)[metadata_neftel$CellAssignment == "Malignant"]]

malignant_metadata_neftel <- metadata_neftel[metadata_neftel$CellAssignment == "Malignant",]

malignant_data_sparse_neftel <- as(as.matrix(malignant_data_neftel), "sparseMatrix")

#################################
#### Step 2: Data processing ####
#################################

# Create CellRouter object for Neftel et al dataset.
cellrouter_neftel <- CreateCellRouter(
    malignant_data_sparse_neftel,
    min.cells = 0,
    min.genes = 0,
    is.expr = 0
)
cellrouter_neftel@assays$RNA@rawdata <- cellrouter_neftel@assays$RNA@rawdata[rownames(cellrouter_neftel@assays$RNA@ndata), colnames(cellrouter_neftel@assays$RNA@ndata)]

# Ensure row names in metadata table have the same order as row names in the sampTab slot of the RNA assay & add
# sample identification to the CellRouter object.
malignant_metadata_neftel <- malignant_metadata_neftel[rownames(cellrouter_neftel@assays$RNA@sampTab),]
cellrouter_neftel@assays$RNA@sampTab$samples <- malignant_metadata_neftel$Sample

# Adding normalized data to ndata slot in the RNA assay.
cellrouter_neftel@assays$RNA@ndata <- malignant_data_sparse_neftel

# Transpose normalized data & ensure row names have the same order as row names in the sampTab slot of the RNA assay.
transposed_norm_data <- as.data.frame(t(as.matrix(cellrouter_neftel@assays$RNA@ndata)))
transposed_norm_data <- transposed_norm_data[rownames(cellrouter_neftel@assays$RNA@sampTab),]

# Concatenate dataframes to have gene expression and metadata information in the same table.
expr <- cbind(t(as.matrix(cellrouter_neftel@assays$RNA@ndata)), cellrouter_neftel@assays$RNA@sampTab)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
samples.neftel <- ggplot(
    expr,
    aes(x = reorder(Sample, PRNP), y = PRNP, fill = Sample)
) +
    geom_boxplot(color = "antiquewhite4", outlier.colour = "antiquewhite4") +
    xlab("Sample") +
    ylab("PRNP expression") +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5)
    )

# Scale data according to all genes.
cellrouter_neftel <- scaleData(cellrouter_neftel, blocksize = nrow(cellrouter_neftel@assays$RNA@ndata))

# Perform PCA.
cellrouter_neftel <- computePCA(cellrouter_neftel, num.pcs = 50, seed = 42)

# Compute UMAP.
cellrouter_neftel <- computeUMAP(cellrouter_neftel, num.pcs = 15, seed = 42)

# Save processed CellRouter object.
save(cellrouter_neftel, file = "results/neftel/Neftel_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")

###############################################################
#### Step 3: Establishing groups with distinct PRNP levels ####
###############################################################

# Getting cells that express PRNP.
prnp_pos_neftel <- colnames(cellrouter_neftel@assays$RNA@ndata)[cellrouter_neftel@assays$RNA@ndata["PRNP",] > 0]
# Getting cells that do not express PRNP.
prnp_neg_neftel <- colnames(cellrouter_neftel@assays$RNA@ndata)[cellrouter_neftel@assays$RNA@ndata["PRNP",] == 0]

cellrouter_neftel@assays$RNA@sampTab <- cellrouter_neftel@assays$RNA@sampTab %>% mutate(PRNP_status = case_when(
    rownames(cellrouter_neftel@assays$RNA@sampTab) %in% prnp_pos_neftel ~ "PRNP_positive_cells",
    rownames(cellrouter_neftel@assays$RNA@sampTab) %in% prnp_neg_neftel ~ "PRNP_negative_cells")
)

# Ensure row names in cell embeddings of UMAP slot have the same order as the row names of the transposed, normalized count data.
cellrouter_neftel@umap$cell.embeddings <- cellrouter_neftel@umap$cell.embeddings[rownames(transposed_norm_data),]

# Concatenate dataframes to have cell embeddings and gene expression information in the same table.
dim.red.neftel <- cbind(cellrouter_neftel@umap$cell.embeddings, transposed_norm_data)
colnames(dim.red.neftel)[1] <- "UMAP1"
colnames(dim.red.neftel)[2] <- "UMAP2"

# Ensure row names have the same order as the row names in the sampTab slot of the RNA assay.
dim.red.neftel <- dim.red.neftel[rownames(cellrouter_neftel@assays$RNA@sampTab),]

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.neftel <- cbind(dim.red.neftel, cellrouter_neftel@assays$RNA@sampTab)

#### @ FIGURE 2B, CENTRAL PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
umap.neftel.prnp <- ggplot(
    dim.red.neftel %>% dplyr::arrange(PRNP),
    aes(x = UMAP1, y = UMAP2, color = PRNP)
) +
    geom_point(size = 0.1) +
    theme(
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key.width = unit(1.2, 'cm'),
        legend.key.height = unit(0.5, 'cm')
    ) +
    scale_color_gradientn(
        colors = colorRampPalette((brewer.pal(9, "RdPu")))(1000),
        limits = c(min(dim.red.neftel$PRNP), max(dim.red.neftel$PRNP))
    )

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample ID.
umap.neftel.sample <- ggplot(
    dim.red.neftel,
    aes(x = UMAP1, y = UMAP2, color = Sample)
) +
    geom_point(size = 0.1) +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key.size = unit(2, 'cm')
    ) + 
    guides(
        color = guide_legend(override.aes = list(size = 0.8),
        keywidth = 0,
        keyheight = 0,
        ncol = 3,
        title = "Sample")
    )

#############################################################
#### Step 4: Finding signatures of PRNP+ and PRNP- cells ####
#############################################################

# Find PRNP+ and PRNP- marker gene signatures & save results to output file.
markers_neftel <- findSignatures(
    object = cellrouter_neftel,
    assay.type = "RNA",
    column = "PRNP_status",
    test.use = "wilcox"
)
write.csv(markers_neftel, "results/neftel/neftel_PRNP+_vs_PRNP-_signature_markers.csv")


# This is a volcano plot that, at the moment, is not shown in the paper.
ggplot(
    markers_neftel,
    aes(x = fc, y = -log10(pval))
) +
    geom_point(size = 1, fill = "royalblue", shape = 21) +
    theme_classic() + 
    theme(
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)
    ) +
    scale_x_break(c(1.8,4.5), ticklabels = c(4.5)) +
    ylim(c(0,100)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ylab(expression(-log[10](p-value))) +
    xlab(expression(log[2]("Fold change")))

######################################
#### Step 5: Correlation analysis ####
######################################

df_neftel <- data.frame(matrix(0, ncol=4))
colnames(df_neftel) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_neftel@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_neftel@assays$RNA@ndata[g,]), as.numeric(cellrouter_neftel@assays$RNA@ndata['PRNP',]))
  df_neftel[i,'gene1'] <- 'PRNP'
  df_neftel[i,'gene2'] <- g
  df_neftel[i,'correlation'] <- c$estimate # Store the correlation.
  df_neftel[i,'pvalue'] <- c$p.value # Store the p-value.
  i <- i + 1
}
write.csv(df_neftel, file = "results/neftel/Neftel_Correlation_PRNPvsAllGenes.csv")

df_neftel <- df_neftel[!is.na(df_neftel$correlation),]
df_neftel <- df_neftel[order(df_neftel$correlation, decreasing = TRUE),]
df_neftel <- df_neftel[-c(1),]
rownames(df_neftel) <- 1:length(rownames(df_neftel))

df_neftel <- df_neftel %>% mutate(significant = case_when(
    pvalue <= 0.05 ~ "Significant",
    pvalue > 0.05 ~ "Non-significant"
))

#### @ FIGURE 2C, CENTRAL PANEL (MAIN) @ #### 
neftel.cor <- ggplot(
    df_neftel,
    aes(y = correlation, x = 1:length(rownames(df_neftel)), color = factor(significant))
) +
    geom_point() +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ylab("Pearson correlation") +
    xlab("Number of genes") +
    scale_color_manual(values = c("gray", "deeppink"), name = "Significance")  +
    ggtitle("Neftel et al") +
    theme(plot.title = element_text(hjust = 0.5))

###################################
###################################
#---- Richards et al. dataset ----#
###################################
###################################

##########################################################
#### Step 1: Loading and subsetting data for analysis ####
##########################################################

GBM_44k_raw_data <- read.csv(
    "data/richards/Richards_NatureCancer_GBM_scRNAseq_counts.csv",
    header = TRUE,
    row.names = 1,
    sep = ",",
    check.names = FALSE
)

GSC_and_whole_tumor_metadata <- read.table(
    "data/richards/GSCs_Tumour_MetaData.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

rownames(GSC_and_whole_tumor_metadata) <- gsub('-', '.', rownames(GSC_and_whole_tumor_metadata))
rownames(GSC_and_whole_tumor_metadata) <- gsub('GBM_', '', rownames(GSC_and_whole_tumor_metadata))

malignant_data_richards <- GBM_44k_raw_data[,colnames(GBM_44k_raw_data) %in% rownames(GSC_and_whole_tumor_metadata)[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR"]]

malignant_metadata_richards <- GSC_and_whole_tumor_metadata[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR",]

malignant_data_sparse_richards <- as(as.matrix(malignant_data_richards), "sparseMatrix")

#################################
#### Step 2: Data processing ####
#################################

# Create CellRouter object for Richards et al dataset.
cellrouter_richards <- CreateCellRouter(
    malignant_data_sparse_richards, 
    min.genes = 0,
    min.cells = 0,
    is.expr = 0
)
cellrouter_richards@assays$RNA@rawdata <- cellrouter_richards@assays$RNA@rawdata[rownames(cellrouter_richards@assays$RNA@ndata), colnames(cellrouter_richards@assays$RNA@ndata)]

# Ensure row names in metadata table have the same order as row names in the sampTab slot of the RNA assay & add
# sample identification to the CellRouter object.
malignant_metadata_richards <- malignant_metadata_richards[rownames(cellrouter_richards@assays$RNA@sampTab),]
cellrouter_richards@assays$RNA@sampTab$samples <- malignant_metadata_richards$Sample

# Normalize malignant count data.
cellrouter_richards <- Normalize(cellrouter_richards)

# Transpose normalized data & ensure row names have the same order as row names in the sampTab slot of the RNA assay.
transposed_norm_data <- as.data.frame(t(as.matrix(cellrouter_richards@assays$RNA@ndata)))
transposed_norm_data <- transposed_norm_data[rownames(cellrouter_richards@assays$RNA@sampTab),]

# Concatenate dataframes to have gene expression and metadata information in the same table.
expr <- cbind(t(as.matrix(cellrouter_richards@assays$RNA@ndata)), cellrouter_richards@assays$RNA@sampTab)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
samples.richards <- ggplot(
    expr,
    aes(x = reorder(Sample, PRNP), y = PRNP, fill = Sample)
) +
    geom_boxplot(color = "antiquewhite4", outlier.colour = "antiquewhite4") +
    xlab("Sample") +
    ylab("PRNP expression") +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5)
    )

# Scale data according to all genes.
cellrouter_richards <- scaleData(cellrouter_richards, blocksize = nrow(cellrouter_neftel@assays$RNA@ndata))

# Perform PCA.
cellrouter_richards <- computePCA(cellrouter_richards, num.pcs = 50, seed = 42)

# Compute UMAP.
cellrouter_richards <- computeUMAP(cellrouter_richards, num.pcs = 15, seed = 42)

# Save processed CellRouter object.
save(cellrouter_richards, file = "results/richards/Richards_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")

###############################################################
#### Step 3: Establishing groups with distinct PRNP levels ####
###############################################################

# Getting cells that express PRNP.
prnp_pos_richards <- colnames(cellrouter_richards@assays$RNA@ndata)[cellrouter_richards@assays$RNA@ndata["PRNP",] > 0]
# Getting cells that do not express PRNP.
prnp_neg_richards <- colnames(cellrouter_richards@assays$RNA@ndata)[cellrouter_richards@assays$RNA@ndata["PRNP",] == 0]

cellrouter_richards@assays$RNA@sampTab <- cellrouter_richards@assays$RNA@sampTab %>% mutate(PRNP_status = case_when(
    rownames(cellrouter_richards@assays$RNA@sampTab) %in% prnp_pos_richards ~ "PRNP_positive_cells",
    rownames(cellrouter_richards@assays$RNA@sampTab) %in% prnp_neg_richards ~ "PRNP_negative_cells")
)

# Ensure row names in cell embeddings of UMAP slot have the same order as the row names of the transposed, normalized count data.
cellrouter_richards@umap$cell.embeddings <- cellrouter_richards@umap$cell.embeddings[rownames(transposed_norm_data),]

# Concatenate dataframes to have cell embeddings and gene expression information in the same table.
dim.red.richards <- cbind(cellrouter_richards@umap$cell.embeddings, transposed_norm_data)
colnames(dim.red.richards)[1] <- "UMAP1"
colnames(dim.red.richards)[2] <- "UMAP2"

# Ensure row names have the same order as the row names in the sampTab slot of the RNA assay.
dim.red.richards <- dim.red.richards[rownames(cellrouter_richards@assays$RNA@sampTab),]

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.richards <- cbind(dim.red.richards, cellrouter_richards@assays$RNA@sampTab)

#### @ FIGURE 2B, RIGHT PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
umap.richards.prnp <- ggplot(
    dim.red.richards %>% dplyr::arrange(PRNP),
    aes(x = UMAP1, y = UMAP2, color = PRNP)
) +
    geom_point(size = 0.1) +
    theme(
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key.width = unit(1.2, 'cm'),
        legend.key.height = unit(0.5, 'cm')
    ) +
    scale_color_gradientn(
        colors = colorRampPalette((brewer.pal(9, "RdPu")))(1000),
        limits = c(min(dim.red.richards$PRNP), max(dim.red.richards$PRNP))
    )

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample ID.
umap.richards.sample <- ggplot(
    dim.red.richards,
    aes(x = UMAP1, y = UMAP2, color = Sample)
) +
    geom_point(size = 0.1) +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key.size = unit(2, 'cm')
    ) + 
    guides(
        color = guide_legend(override.aes = list(size = 0.8),
        keywidth = 0,
        keyheight = 0,
        ncol = 3,
        title = "Sample")
    )
 
#############################################################
#### Step 4: Finding signatures of PRNP+ and PRNP- cells ####
#############################################################

# Find PRNP+ and PRNP- marker gene signatures & save results to output file.
markers_richards <- findSignatures(
    object = cellrouter_richards, 
    assay.type = "RNA",
    column = "PRNP_status", 
    test.use = "wilcox"
)
write.csv(markers_richards, "results/richards/Richards_PRNP+_vs_PRNP-_signature_markers.csv")

# This is a volcano plot that, at the moment, is not shown in the paper.
ggplot(
    markers_richards,
    aes(x = fc, y = -log10(pval))
) +
    geom_point(size = 1, fill = "royalblue", shape = 21) +
    theme_classic() + 
    theme(
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)
        ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ylab(expression(-log[10](p-value))) +
    xlab(expression(log[2]("Fold change")))

######################################
#### Step 5: Correlation analysis ####
######################################

df_richards <- data.frame(matrix(0, ncol=4))
colnames(df_richards) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_richards@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_richards@assays$RNA@ndata[g,]), as.numeric(cellrouter_richards@assays$RNA@ndata['PRNP',]))
  df_richards[i,'gene1'] <- 'PRNP'
  df_richards[i,'gene2'] <- g
  df_richards[i,'correlation'] <- c$estimate # Store the correlation.
  df_richards[i,'pvalue'] <- c$p.value # Store the p-value.
  i <- i + 1
}
write.csv(df_richards, file = "results/richards/Richards_Correlation_PRNPvsAllGenes.csv")

df_richards <- df_richards[!is.na(df_richards$correlation),]
df_richards <- df_richards[order(df_richards$correlation, decreasing = TRUE),]
df_richards <- df_richards[-c(1),]
rownames(df_richards) <- 1:length(rownames(df_richards))

df_richards <- df_richards %>% mutate(significant = case_when(
    pvalue <= 0.05 ~ "Significant",
    pvalue > 0.05 ~ "Non-significant"
))

#### @ FIGURE 2C, RIGHT PANEL (MAIN) @ #### 
richards.cor <- ggplot(
    df_richards,
    aes(y = correlation, x = 1:length(rownames(df_richards)), color = factor(significant))
) +
    geom_point() +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ylab("Pearson correlation") +
    xlab("Number of genes") +
    scale_color_manual(values = c("gray", "deeppink"), name = "Significance")  +
    ggtitle("Richards et al") +
    theme(plot.title = element_text(hjust = 0.5))

######################################

sessionInfo()