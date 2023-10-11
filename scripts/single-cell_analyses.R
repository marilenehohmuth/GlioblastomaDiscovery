# SINGLE-CELL ANALYSES FOR PAPER ###############################################

# R version: 4.0.3
# VERSIONS
library(fusca)                    # 1.2.1
library(ggplot2)                  # 3.3.5
library(dplyr)                    # 1.0.7
library(Matrix)                   # 1.3-4
library(RColorBrewer)             # 1.1-2
library(Hmisc)                    # 4.5-0
library(corrplot)                 # 0.92

setwd("manuscript2022/")

# PART 1 - ANALYSIS OF GBM DATA FROM DARMANIS ET AL DATASET ####################

########## Step 1: Loading and subsetting data for analysis ##########

# Load Darmanis et al count data.
data_darmanis <- read.csv("data/darmanis/GBM_raw_gene_counts.csv",
                          header = TRUE,
                          row.names = 1,
                          sep = " ",
                          check.names = FALSE)

# Load Darmanis et al metadata.
metadata_darmanis <- read.csv("data/darmanis/GBM_metadata.csv",
                              header = TRUE,
                              row.names = 1,
                              sep = " ")

# Subset count data to keep only malignant cells.
malignant_data_darmanis <- data_darmanis[,colnames(data_darmanis) %in%
                                           rownames(metadata_darmanis)[metadata_darmanis$Cluster_2d %in% c(1, 4, 11)]]

# Subset metadata to keep only malignant cells.
malignant_metadata_darmanis <- metadata_darmanis[metadata_darmanis$Cluster_2d %in% c(1, 4, 11),]

# Convert count data data.frame to sparse matrix.
malignant_data_sparse_darmanis <- as(as.matrix(malignant_data_darmanis), "sparseMatrix")

########## Step 2: Data processing ##########

# Create CellRouter object for Darmanis et al dataset.
cellrouter_darmanis <- CreateCellRouter(malignant_data_sparse_darmanis,
                                        min.genes = 0,
                                        min.cells = 0,
                                        is.expr = 0)

cellrouter_darmanis@assays$RNA@rawdata <- cellrouter_darmanis@assays$RNA@rawdata[rownames(cellrouter_darmanis@assays$RNA@ndata), 
                                                                                 colnames(cellrouter_darmanis@assays$RNA@ndata)]

# Add sample identification to CellRouter object.
cellrouter_darmanis@assays$RNA@sampTab$samples <- malignant_metadata_darmanis$Sample.name
# Note: row order of malignant_metadata_darmanis = row order of sampTab in the RNA assay.

# Normalize malignant count data.
cellrouter_darmanis <- Normalize(cellrouter_darmanis)

expr <- cbind(t(as.matrix(cellrouter_darmanis@assays$RNA@ndata)), 
              cellrouter_darmanis@assays$RNA@sampTab)
# Note: row order of of sampTab in the RNA assay = row order of transposed ndata.

boxplot.samples.darmanis <- ggplot(expr,
                                   aes(x = reorder(samples, PRNP),
                                   y = PRNP,
                                   fill = samples)) +
                              geom_boxplot(color = "antiquewhite4", 
                                           outlier.colour = "antiquewhite4") +
                              ylab("Normalized PRNP counts") +
                              xlab("Sample") +
                              theme_classic() +
                              theme(legend.position = "none",
                                    axis.text.x = element_text(size = 5, 
                                                               angle = 45,
                                                               vjust = 1,
                                                               hjust = 1))

# Scale data according to all genes.
cellrouter_darmanis <- scaleData(cellrouter_darmanis,
                                 blocksize = nrow(cellrouter_darmanis@assays$RNA@ndata))

# Perform PCA.
cellrouter_darmanis <- computePCA(cellrouter_darmanis, 
                                  num.pcs = 50, 
                                  seed = 42)

# Compute UMAP.
cellrouter_darmanis <- computeUMAP(cellrouter_darmanis, 
                                   num.pcs = 15, 
                                   seed = 42)

# Save processed CellRouter object of the Darmanis et al dataset.
save(cellrouter_darmanis,
     file = "results/darmanis/Darmanis_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")

########## Step 3: Establishing groups with distinct PRNP levels ##########

prnp_pos_darmanis <- colnames(cellrouter_darmanis@assays$RNA@ndata)[cellrouter_darmanis@assays$RNA@ndata["PRNP",] > 0]
prnp_neg_darmanis <- colnames(cellrouter_darmanis@assays$RNA@ndata)[cellrouter_darmanis@assays$RNA@ndata["PRNP",] == 0]

cellrouter_darmanis@assays$RNA@sampTab <- cellrouter_darmanis@assays$RNA@sampTab %>% 
  mutate(PRNP_status = case_when(rownames(cellrouter_darmanis@assays$RNA@sampTab) %in% prnp_pos_darmanis ~"PRNP_positive_cells",
                                 rownames(cellrouter_darmanis@assays$RNA@sampTab) %in% prnp_neg_darmanis ~"PRNP_negative_cells"))
print("Classified malignant GBM cells according to PRNP expression.")

dim.red.darmanis <- cbind(cellrouter_darmanis@umap$cell.embeddings, 
                          as.data.frame(t(as.matrix(cellrouter_darmanis@assays$RNA@ndata))))
colnames(dim.red.darmanis)[1] <- "UMAP1"
colnames(dim.red.darmanis)[2] <- "UMAP2"
# Note: row order of embeddings = row order of transposed ndata.

dim.red.darmanis <- cbind(dim.red.darmanis,
                          cellrouter_darmanis@assays$RNA@sampTab)
# Note: row order of dim.red.darmanis = row order of sampTab.

ggplot(dim.red.darmanis %>% 
         dplyr::arrange(PRNP),
       aes(x = UMAP1,
           y = UMAP2,
           color = PRNP)) +
  geom_point(size = 1) +
  theme_classic() +
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        #legend.key.width = unit(0.6, 'cm'),
        #legend.key.height = unit(1.988, 'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(1000),
                        limits = c(min(dim.red.darmanis$PRNP),
                                   max(dim.red.darmanis$PRNP)))

umap_sample <- ggplot(dim.red.darmanis,
                      aes(x = UMAP1,
                      y = UMAP2,
                      color = samples)) +
  geom_point(size = 1) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.size = unit(2, 'cm'),
        legend.box.just = "center") + 
  guides(color = guide_legend(override.aes = list(size=3),
         keywidth = 0,
         keyheight = 0,
         ncol = 4,
         title = "")) 

########## Step 4: Finding signatures of PRNP+ and PRNP- cells ##########

# Find PRNP+ and PRNP- marker gene signatures.
markers_darmanis <- findSignatures(object = cellrouter_darmanis,
                                   assay.type = "RNA",
                                   column = "PRNP_status",
                                   test.use = "wilcox")

write.csv(markers_darmanis,
          "results/darmanis/Darmanis_PRNP+_vs_PRNP-_signature_markers.csv")

ggplot(markers_darmanis,
       aes(x = fc,
           y = -log10(pval))) +
  geom_point(size = 1,
             fill = "royalblue",
             shape = 21) +
  theme_classic() + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_y_break(c(50,185), 
                ticklabels=c(190)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ylab(expression(-log[10](p-value))) +
  xlab(expression(log[2]("Fold change")))

########## Step 5: Correlation analyses ##########

print("Starting correlation analysis.")
df_darmanis <- data.frame(matrix(0, ncol=4))
colnames(df_darmanis) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_darmanis@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_darmanis@assays$RNA@ndata[g,]), 
                as.numeric(cellrouter_darmanis@assays$RNA@ndata['PRNP',]))
  df_darmanis[i,'gene1'] <- 'PRNP'
  df_darmanis[i,'gene2'] <- g
  df_darmanis[i,'correlation'] <- c$estimate #stores the correlation
  df_darmanis[i,'pvalue'] <- c$p.value #stores the correlation
  i <- i + 1
}

print("Calculated correlation between PRNP and all genes.")

write.csv(df_darmanis, 
          file = "results/darmanis/Darmanis_Correlation_PRNPvsAllGenes.csv")

# PART 2 - ANALYSIS OF GBM DATA FROM NEFTEL ET AL DATASET ####################

########## Step 1: Loading and subsetting data for analysis ##########

# Load Neftel et al count data.
data_neftel <- read.table("data/neftel/IDHwtGBM.processed.SS2.logTPM.txt",
                        header = TRUE,
                        row.names = 1,
                        sep = "\t",
                        check.names = FALSE)

# Load Neftel et al metadata.
metadata_neftel <- read.table("data/neftel/IDHwt.GBM.Metadata.SS2.txt",
                              header = TRUE,
                              row.names = 1,
                              sep = "\t",
                              check.names = FALSE)

# Subset count data to keep only malignant cells.
malignant_data_neftel <- data_neftel[,colnames(data_neftel) %in% 
                                       rownames(metadata_neftel)[metadata_neftel$CellAssignment == "Malignant"]]

malignant_metadata_neftel <- metadata_neftel[metadata_neftel$CellAssignment == "Malignant",]

malignant_data_sparse_neftel <- as(as.matrix(malignant_data_neftel), "sparseMatrix")

########## Step 2: Data processing ##########

# Create CellRouter object for Neftel et al dataset.
cellrouter_neftel <- CreateCellRouter(malignant_data_sparse_neftel,
                                      min.cells = 0,
                                      min.genes = 0,
                                      is.expr = 0)
print("Created CellRouter object for Neftel dataset.")

cellrouter_neftel@assays$RNA@rawdata <- cellrouter_neftel@assays$RNA@rawdata[rownames(cellrouter_neftel@assays$RNA@ndata), 
                                                                             colnames(cellrouter_neftel@assays$RNA@ndata)]

# Reorder malignant metadata to match row order of the sampTab in the RNA assay.
malignant_metadata_neftel <- malignant_metadata_neftel[rownames(cellrouter_neftel@assays$RNA@sampTab),]

cellrouter_neftel@assays$RNA@sampTab$samples <- malignant_metadata_neftel$Sample

cellrouter_neftel@assays$RNA@ndata <- malignant_data_sparse_neftel

expr <- cbind(t(as.matrix(cellrouter_neftel@assays$RNA@ndata)), 
              cellrouter_neftel@assays$RNA@sampTab)
# Note: row order of transposed ndata = row order of sampTab.

samples.neftel <- ggplot(expr,
                         aes(x = reorder(Sample, PRNP),
                             y = PRNP,
                             fill = Sample)) +
                  geom_boxplot(color = "antiquewhite4", 
                               outlier.colour = "antiquewhite4") +
                  ylab("PRNP expression") +
                  xlab("Sample") +
                  theme_classic() +
                  theme(legend.position = "none",
                        axis.text.x = element_text(angle = 45,
                                                   hjust = 1,
                                                   vjust = 1,
                                                   size = 5))

# Scale data according to all genes.
cellrouter_neftel <- scaleData(cellrouter_neftel,
                               blocksize = nrow(cellrouter_neftel@assays$RNA@ndata))

# Perform PCA.
cellrouter_neftel <- computePCA(cellrouter_neftel, num.pcs = 50, seed = 42)

# Compute UMAP.
cellrouter_neftel <- computeUMAP(cellrouter_neftel, num.pcs = 15, seed = 42)

# Save processed CellRouter object.
save(cellrouter_neftel,
     file = "results/neftel/Neftel_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")

########## Step 3: Establishing groups with distinct PRNP levels ##########

prnp_pos_neftel <- colnames(cellrouter_neftel@assays$RNA@ndata)[cellrouter_neftel@assays$RNA@ndata["PRNP",] > 0]
prnp_neg_neftel <- colnames(cellrouter_neftel@assays$RNA@ndata)[cellrouter_neftel@assays$RNA@ndata["PRNP",] == 0]

cellrouter_neftel@assays$RNA@sampTab <- cellrouter_neftel@assays$RNA@sampTab %>% 
  mutate(PRNP_status = case_when(rownames(cellrouter_neftel@assays$RNA@sampTab) %in% prnp_pos_neftel ~"PRNP_positive_cells",
                                 rownames(cellrouter_neftel@assays$RNA@sampTab) %in% prnp_neg_neftel ~"PRNP_negative_cells"))


dim.red.neftel <- cbind(cellrouter_neftel@umap$cell.embeddings, 
                        as.data.frame(t(as.matrix(cellrouter_neftel@assays$RNA@ndata))))
colnames(dim.red.neftel)[1] <- "UMAP1"
colnames(dim.red.neftel)[2] <- "UMAP2"
# Note: row order of embeddings = row order of ndata.

dim.red.neftel <- cbind(dim.red.neftel,
                        cellrouter_neftel@assays$RNA@sampTab)


ggplot(dim.red.neftel %>% 
         dplyr::arrange(PRNP),
       aes(x = UMAP1,
           y = UMAP2,
           color = PRNP)) +
  geom_point(size = 0.1) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key.width = unit(1.2, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_gradientn(colors = colorRampPalette((brewer.pal(9, "RdPu")))(1000),
                        limits = c(min(dim.red.neftel$PRNP),
                                   max(dim.red.neftel$PRNP)))

umap.neftel.sample <- ggplot(dim.red.neftel,
                               aes(x = UMAP1,
                                   y = UMAP2,
                                   color = Sample)) +
                      geom_point(size = 0.1) +
                      theme(legend.position = "bottom",
                            legend.text = element_text(size=5),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            legend.key.size = unit(2, 'cm')) + 
                      guides(color = guide_legend(override.aes = list(size=0.8),
                                                  keywidth = 0,
                                                  keyheight = 0,
                                                  ncol = 3,
                                                  title = "Sample"))

print("Plotted UMAP displaying PRNP expression across cells.")

########## Step 4: Finding signatures of PRNP+ and PRNP- cells ##########

markers_neftel <- findSignatures(object = cellrouter_neftel,
                                 assay.type = "RNA",
                                 column = "PRNP_status",
                                 test.use = "wilcox")
print("Found PRNP+ and PRNP- marker gene signatures.")

write.csv(markers_neftel,
          "results/neftel/neftel_PRNP+_vs_PRNP-_signature_markers.csv")

ggplot(markers_neftel,
       aes(x = fc,
           y = -log10(pval))) +
  geom_point(size = 1,
             fill = "royalblue",
             shape = 21) +
  theme_classic() + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_x_break(c(1.8,4.5), 
                ticklabels=c(4.5)) +
  ylim(c(0,100)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ylab(expression(-log[10](p-value))) +
  xlab(expression(log[2]("Fold change")))

########## Step 5: Correlation analyses ##########

# Calculating the Pearson correlation between PRNP expression and the expression 
# levels of each gene expressed by malignant cells in the Neftel et al dataset.
df_neftel <- data.frame(matrix(0, ncol=4))
colnames(df_neftel) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_neftel@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_neftel@assays$RNA@ndata[g,]), 
                as.numeric(cellrouter_neftel@assays$RNA@ndata['PRNP',]))
  df_neftel[i,'gene1'] <- 'PRNP'
  df_neftel[i,'gene2'] <- g
  df_neftel[i,'correlation'] <- c$estimate #stores the correlation
  df_neftel[i,'pvalue'] <- c$p.value #stores the correlation
  i <- i + 1
}

write.csv(df_neftel, 
          file = "results/neftel/Neftel_Correlation_PRNPvsAllGenes.csv")

# PART 4 - ANALYSIS OF GBM DATA FROM RICHARDS ET AL DATASET ####################

########## Step 1: Loading and subsetting data for analysis ##########

GBM_44k_raw_data <- read.csv("data/richards/Richards_NatureCancer_GBM_scRNAseq_counts.csv",
                             header = TRUE,
                             row.names = 1,
                             sep = ",",
                             check.names = FALSE)
print("Loaded Richards count data.")

GSC_and_whole_tumor_metadata <- read.table("data/richards/GSCs_Tumour_MetaData.txt",
                                      header = TRUE,
                                      row.names = 1,
                                      sep = "\t",
                                      check.names = FALSE)
print("Loaded Richards metadata.")

rownames(GSC_and_whole_tumor_metadata) <- gsub('-', '.', rownames(GSC_and_whole_tumor_metadata))
rownames(GSC_and_whole_tumor_metadata) <- gsub('GBM_', '', rownames(GSC_and_whole_tumor_metadata))

malignant_data_richards <- GBM_44k_raw_data[,colnames(GBM_44k_raw_data) %in%
                                              rownames(GSC_and_whole_tumor_metadata)[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR"]]

malignant_metadata_richards <- GSC_and_whole_tumor_metadata[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR",]

print("Filtered raw data.")

malignant_data_sparse_richards <- as(as.matrix(malignant_data_richards), "sparseMatrix")

########## Step 2: Data processing ##########

# Create CellRouter object for Richards et al dataset.
cellrouter_richards <- CreateCellRouter(malignant_data_sparse_richards, 
                                        min.genes = 0,
                                        min.cells = 0,
                                        is.expr = 0)

cellrouter_richards@assays$RNA@rawdata <- cellrouter_richards@assays$RNA@rawdata[rownames(cellrouter_richards@assays$RNA@ndata), 
                                                                                 colnames(cellrouter_richards@assays$RNA@ndata)]

cellrouter_richards@assays$RNA@sampTab$samples <- malignant_metadata_richards$Sample.ID
# Note: row order of sampTab = row order of metadata.

# Normalize data.
cellrouter_richards <- Normalize(cellrouter_richards)

expr <- cbind(t(as.matrix(cellrouter_richards@assays$RNA@ndata)), 
              cellrouter_richards@assays$RNA@sampTab)
# Note: row order of transposed ndata = row order of sampTab.

samples.richards <- ggplot(expr,
                           aes(x = reorder(Sample.ID, PRNP),
                           y = PRNP,
                           fill = Sample.ID)) +
                    geom_boxplot(color = "antiquewhite4", 
                                outlier.colour = "antiquewhite4") +
                    ylab("PRNP expression") +
                    xlab("Sample") +
                    theme_classic() +
                    theme(legend.position = "none",
                          axis.text.x = element_text(angle = 45,
                                                     hjust = 1,
                                                     vjust = 1,
                                                     size = 5))

cellrouter_richards <- scaleData(cellrouter_richards, 
                                 blocksize = nrow(cellrouter_richards@assays$RNA@ndata))
print("Scaled data according to all genes.")

cellrouter_richards <- computePCA(cellrouter_richards, num.pcs = 50, seed = 42)
print("Performed PCA.")

cellrouter_richards <- computeUMAP(cellrouter_richards, num.pcs = 15, seed = 42)
print("Computed UMAP.")

save(cellrouter_richards, 
     file = "results/richards/Richards_Normalized_Scaled_PCA_UMAP_CellRouter_object.R")
print("Saved processed CellRouter object.")

########## Step 3: Establishing groups with distinct PRNP levels ##########

prnp_pos_richards <- colnames(cellrouter_richards@assays$RNA@ndata)[cellrouter_richards@assays$RNA@ndata["PRNP",] > 0]
prnp_neg_richards <- colnames(cellrouter_richards@assays$RNA@ndata)[cellrouter_richards@assays$RNA@ndata["PRNP",] == 0]

cellrouter_richards@assays$RNA@sampTab <- cellrouter_richards@assays$RNA@sampTab %>% 
  mutate(PRNP_status = case_when(rownames(cellrouter_richards@assays$RNA@sampTab) %in% prnp_pos_richards ~"PRNP_positive_cells",
                                 rownames(cellrouter_richards@assays$RNA@sampTab) %in% prnp_neg_richards ~"PRNP_negative_cells"))

dim.red.richards <- cbind(cellrouter_richards@umap$cell.embeddings, 
                        as.data.frame(t(as.matrix(cellrouter_richards@assays$RNA@ndata))))
colnames(dim.red.richards)[1] <- "UMAP1"
colnames(dim.red.richards)[2] <- "UMAP2"
# identical rows

dim.red.richards <- cbind(dim.red.richards,
                          cellrouter_richards@assays$RNA@sampTab)

umap.richards.prnp <- ggplot(dim.red.richards %>% 
                               dplyr::arrange(PRNP),
                             aes(x = UMAP1,
                                 y = UMAP2,
                                 color = PRNP)) +
                      geom_point(size = 0.1) +
                      theme(legend.position = "bottom",
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            legend.key.width = unit(1.2, 'cm'),
                            legend.key.height = unit(0.5, 'cm')) +
                      scale_color_gradientn(colors = colorRampPalette((brewer.pal(9, "RdPu")))(1000),
                                            limits = c(min(dim.red.richards$PRNP),
                                                       max(dim.red.richards$PRNP)))


umap.richards.sample <- ggplot(dim.red.richards,
                               aes(x = UMAP1,
                                   y = UMAP2,
                                   color = Sample.ID)) +
                        geom_point(size = 0.1) +
                        theme(legend.position = "bottom",
                              legend.text = element_text(size=5),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                              legend.key.size = unit(2, 'cm')) + 
                        guides(color = guide_legend(override.aes = list(size=0.8),
                               keywidth = 0,
                               keyheight = 0,
                               ncol = 3,
                               title = "Sample",
                               title.position = "left"))
 
########## Step 4: Finding signatures of PRNP+ and PRNP- cells ##########

markers_Richards <- findSignatures(object = cellrouter_richards, 
                                   assay.type = "RNA",
                                   column = "PRNP_status", 
                                   test.use = "wilcox")
print("Found PRNP+ and PRNP- marker gene signatures.")

write.csv(markers_richards, 
          "results/richards/Richards_PRNP+_vs_PRNP-_signature_markers.csv")

ggplot(markers_richards,
       aes(x = fc,
           y = -log10(pval))) +
  geom_point(size = 1,
             fill = "royalblue",
             shape = 21) +
  theme_classic() + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  # scale_y_break(c(50,185), 
  #               ticklabels=c(190)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ylab(expression(-log[10](p-value))) +
  xlab(expression(log[2]("Fold change")))

########## Step 5: Correlation analyses ##########

print("Starting correlation analysis.")
df_richards <- data.frame(matrix(0, ncol=4))
colnames(df_richards) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i=1
for(g in rownames(cellrouter_richards@assays$RNA@ndata)){
  c <- cor.test(as.numeric(cellrouter_richards@assays$RNA@ndata[g,]), 
                as.numeric(cellrouter_richards@assays$RNA@ndata['PRNP',]))
  df_richards[i,'gene1'] <- 'PRNP'
  df_richards[i,'gene2'] <- g
  df_richards[i,'correlation'] <- c$estimate #stores the correlation
  df_richards[i,'pvalue'] <- c$p.value #stores the correlation
  i <- i + 1
}

write.csv(df_richards, 
          file = "results/richards/richards_Correlation_PRNPvsAllGenes.csv")

# PART 5 - DARMANIS + NEFTEL + RICHARDS ####################


df_darmanis <- subset(df_darmanis,
                      is.na(df_darmanis$correlation) == FALSE)

df_darmanis <- df_darmanis[order(df_darmanis$correlation, decreasing = TRUE),]
df_darmanis <- df_darmanis[-c(1),]
rownames(df_darmanis) <- 1:length(rownames(df_darmanis))

df_darmanis <- df_darmanis %>%
  mutate(significant = case_when(pvalue <= 0.05 ~ "Significant",
                                 pvalue > 0.05 ~ "Non-significant"))

darmanis.cor <- ggplot(df_darmanis,
                       aes(y = correlation,
                           x = 1:length(rownames(df_darmanis)),
                           color = factor(significant))) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0,
             linetype = 'dotted') +
  ylab("Spearman correlation") +
  xlab("Number of genes") +
  scale_color_manual(values = c("gray", "deeppink"),
                     name = "Significance")  +
  ggtitle("Darmanis et al") +
  theme(plot.title = element_text(hjust = 0.5))


df_neftel <- subset(df_neftel,
                    is.na(df_neftel$correlation) == FALSE)

df_neftel <- df_neftel[order(df_neftel$correlation, decreasing = TRUE),]
df_neftel <- df_neftel[-c(1),]
rownames(df_neftel) <- 1:length(rownames(df_neftel))

df_neftel <- df_neftel %>%
  mutate(significant = case_when(pvalue <= 0.05 ~ "Significant",
                                 pvalue > 0.05 ~ "Non-significant"))

neftel.cor <- ggplot(df_neftel,
                     aes(y = correlation,
                         x = 1:length(rownames(df_neftel)),
                         color = factor(significant))) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0,
             linetype = 'dotted') +
  ylab("Spearman correlation") +
  xlab("Number of genes") +
  scale_color_manual(values = c("gray", "deeppink"),
                     name = "Significance")  +
  ggtitle("Neftel et al") +
  theme(plot.title = element_text(hjust = 0.5))

df_richards <- subset(df_richards,
                      is.na(df_richards$correlation) == FALSE)

df_richards <- df_richards[order(df_richards$correlation, decreasing = TRUE),]
df_richards <- df_richards[-c(1),]
rownames(df_richards) <- 1:length(rownames(df_richards))

df_richards <- df_richards %>%
  mutate(significant = case_when(pvalue <= 0.05 ~ "Significant",
                                 pvalue > 0.05 ~ "Non-significant"))

richards.cor <- ggplot(df_richards,
                       aes(y = correlation,
                           x = 1:length(rownames(df_richards)),
                           color = factor(significant))) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0,
             linetype = 'dotted') +
  ylab("Spearman correlation") +
  xlab("Number of genes") +
  scale_color_manual(values = c("gray", "deeppink"),
                     name = "Significance") +
  ggtitle("Richards et al") +
  theme(plot.title = element_text(hjust = 0.5))

corr.intersect <- read.table("results/intersection_correlations/darmanis_significant_correlations_panther_protein_classes.txt",
                             sep = "\t",
                             header = FALSE,
                             row.names = 1)

corr <- ggplot(corr.intersect,
       aes(y = V3,
           x = reorder(V2, -V3),
           fill = V2)) +
  geom_col(width = 1, 
           color = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.position = "none") +
  ylab("Number of proteins") +
  xlab("Protein class") +
  scale_fill_manual(values = qualitative_hcl(length(rownames(corr.intersect)), palette = "pastel1")) 

pdf("results/Darmanis+Neftel+Richards_Correlation_PRNPxEachGene.pdf",
    width = 10,
    height = 3)
ggarrange(plotlist = list(darmanis.cor, neftel.cor, richards.cor),
          align = "hv",
          common.legend = TRUE,
          legend = "bottom",
          ncol = 3,
          nrow = 1)
dev.off()

sessionInfo()


