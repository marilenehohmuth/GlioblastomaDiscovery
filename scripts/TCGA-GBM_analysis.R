# TCGA ANALYSES FOR PAPER ######################################################

# R version: 4.0.3
                                  # VERSIONS
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

print("Loaded packages.")

# setwd("manuscript2022/")

# PART 1 - ANALYSIS OF GBM DATA FROM TCGA ######################################

print("Starting the analysis of GBM data from TCGA.")

########## Step 1: Data retrieval ##########

# Retrieving data from the TCGA-GBM project.
query_GBM <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",
                      experimental.strategy = "RNA-Seq")
                      # legacy = FALSE)
GDCdownload(query_GBM,
            method = "api",
            files.per.chunk = 10)
data_GBM <- GDCprepare(query = query_GBM)

# Raw count data.
GBM_count_data <- as.data.frame(assay(data_GBM))

# Metadata.
GBM_metadata <- as.data.frame(colData(data_GBM))

print("Retrieved GBM data from TCGA.")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

# ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
# ensembl <- useDataset("hsapiens_gene_ensembl",
#                      mart = ensembl)
# genemap <- getBM(filters = c("ensembl_gene_id"),
#                  attributes = c("ensembl_gene_id","hgnc_symbol"),
#                  values = rownames(GBM_count_data),
#                  mart = ensembl)
# idx <- match(rownames(GBM_count_data), genemap$ensembl_gene_id)
# GBM_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]
# 
# rownames(GBM_count_data) <- make.names(GBM_count_data$hgnc_symbol, unique = TRUE)
# GBM_count_data$hgnc_symbol <- NULL
  
########## Step 3: Exploring PRNP expression in different tissue types ##########

# Normalizing raw count data of all samples, tumor and non-tumor (CPM normalization).
tissue_cpm <- cpm(GBM_count_data, 
                  log = FALSE)

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

saveRDS(tissue_log_metadata, 
        file = "results/TCGA/Tissue_logData+Metadata.RDS")

########## Step 4: Selecting primary tumors ##########

# Filtering metadata.
GBM_primary_tumor_metadata <- subset(GBM_metadata, 
                                     definition == "Primary solid Tumor")

# Filtering raw count data.
GBM_primary_tumor_count_data <- GBM_count_data %>% 
  dplyr::select(intersect(rownames(GBM_primary_tumor_metadata), 
                   colnames(GBM_count_data)))

print("Selected primary tumor data.")

########## Step 5: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pGBM_cpm <- cpm(GBM_primary_tumor_count_data, 
                log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pGBM_cpm1 <- pGBM_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pGBM_log <- as.data.frame(log10(pGBM_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: Exploring PRNP levels across different conditions in primary tumors ##########

# Ensuring that the row names of the data.frame with log10(CPM+1) counts of primary tumors follow 
# the same order as the row names of the data.frame with metadata.
pGBM_log.t <- as.data.frame(t(pGBM_log))
pGBM_log.t <- pGBM_log.t[rownames(GBM_primary_tumor_metadata),]

# Mixing data.frames.
pGBM_log_metadata <- cbind(GBM_primary_tumor_metadata, pGBM_log.t)

saveRDS(pGBM_log_metadata, 
        file = "results/TCGA/PrimaryGBMs_logData+Metadata.RDS")


# Plotting PRNP levels across IDH-WT and IDH-mutant primary GBMs.
pGBM.idh.prnp <- ggplot(pGBM_log_metadata,
       aes(x = paper_IDH.status,
           y = ENSG00000171867.17, # PRNP gene
           fill = paper_IDH.status)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ggtitle("Primary GBMs") +
  ylab(expression(paste(italic("PRNP"), " expression"))) +
  scale_x_discrete(labels = c("IDH-mutant", "IDHwt", "Unclassified")) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#E5E0E4")) +
  stat_compare_means(label = "p.format",
                     label.y.npc = "bottom",
                     label.x.npc = "left",
                     size = 4) # Kruskal-Wallis test


pdf("results/TCGA/PRNP_expression_GBM_Mutant_WT.pdf",
    width = 4,
    height = 4.5)
pGBM.idh.prnp
dev.off()

# Plotting PRNP levels across primary GBM subtypes.
## remove NE classification
pGBM_log_metadata <- pGBM_log_metadata %>% 
  mutate(paper_Transcriptome.Subtype_clean = ifelse(paper_Transcriptome.Subtype == 'NE', 
                                                    NA, as.character(paper_Transcriptome.Subtype)))

pGBM.sub.prnp <- ggplot(pGBM_log_metadata,
       aes(x = paper_Transcriptome.Subtype_clean,
           y = ENSG00000171867.17, # PRNP gene
           fill = paper_Transcriptome.Subtype_clean)) + # paper_Transcriptome.Subtype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ylab("") + 
  stat_compare_means(label = "p.format",
                     label.y.npc = "bottom",
                     label.x.npc = "left",
                     size = 4) +
  scale_x_discrete(labels = c("Classical", "Mesenchymal", "Proneural", "Unclassified")) +
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#C77CFF", "#E5E0E4"))
  # scale_x_discrete(labels = c("Classical", "Mesenchymal", "Neural", "Proneural", "Unclassified")) +
  # scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#E5E0E4"))

pdf("results/TCGA/PRNP_expression_across_GBM_subtypes.pdf",
    width = 4,
    height = 4.5)
pGBM.sub.prnp
dev.off()


pdf("results/TCGA/arranged_plots_pGBM.pdf",
    width = 8,
    height = 4)
cowplot::plot_grid(plotlist = list(pGBM.idh.prnp, pGBM.sub.prnp),
                   align = "hv",
                   ncol = 2,
                   nrow = 1,
                   axis = "btlr",
                   rel_widths = c(1.8,2))
dev.off()

########## Step 7: Establishing groups with distinct PRNP levels ##########

PRNP_counts_all <- pGBM_log["ENSG00000171867.17",]
PRNP_counts_all <- as.data.frame(t(PRNP_counts_all))

PRNP_quartiles <- quantile(PRNP_counts_all$ENSG00000171867.17, c(0.25, 0.75))
PRNP_quartiles <- as.data.frame(PRNP_quartiles)
rownames(PRNP_quartiles) <- c("First", "Third")

PRNP_counts_all <- PRNP_counts_all %>% 
  dplyr::mutate(PRNP_status = case_when(ENSG00000171867.17 <= PRNP_quartiles["First",] ~ "PRNP-Low",
                                 ENSG00000171867.17 > PRNP_quartiles["Third",] ~ "PRNP-High"))

PRNP_counts <- PRNP_counts_all[which(PRNP_counts_all$PRNP_status %in% 
                                   c("PRNP-Low","PRNP-High")),]

# This will be the raw count matrix for DESeq2.
PRNP_quartiles_count_data <- GBM_primary_tumor_count_data %>% 
  dplyr::select(rownames(PRNP_counts))

# This will be the metadata for DESeq2.
PRNP_quartiles_status <- PRNP_counts %>% 
  dplyr::select(PRNP_status)

print("Established groups with distinct PRNP levels.")

########## Step 8: Characterizing the PRNP-High and PRNP-Low groups ##########

# Selecting log10(CPM+1) counts of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_log_data <- as.data.frame(pGBM_log) %>% dplyr::select(rownames(PRNP_counts))

# Mixing data.frames.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_status, t(PRNP_quartiles_log_data)) 

# Selecting metadata of primary GBMs classified as PRNP-High and PRNP-Low.
PRNP_quartiles_metadata <- as.data.frame(t(GBM_primary_tumor_metadata)) %>% 
  dplyr::select(rownames(PRNP_counts))

# Mixing data.frames.
PRNP_quartiles_log_metadata <- cbind(PRNP_quartiles_log_metadata, t(PRNP_quartiles_metadata)) 


groups.prnp.idh <- ggplot(PRNP_quartiles_log_metadata,
       aes(x = reorder(PRNP_status, ENSG00000171867.17),
           fill = unlist(paper_IDH.status))) +
  geom_bar(position = "fill",
           color = "black") +
  theme_classic() +
  xlab("Group") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 0, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  xlab("") +
  scale_fill_manual(name = "IDH status",
                    values = c("#00BFC4", "#F8766D", "snow2"),
                    labels = c("IDH-mutant", "IDHwt", "Unclassified")) +
  scale_x_discrete(labels = c(expression(italic("PRNP")^"low"), expression(italic("PRNP")^"high")))

pdf("results/TCGA/IDH_status_across_PRNP-High_and_PRNP-Low_groups.pdf",
    width = 4,
    height = 4.5)
groups.prnp.idh
dev.off()

groups.prnp.sub <- ggplot(PRNP_quartiles_log_metadata,
       aes(x = reorder(PRNP_status, ENSG00000171867.17),
           fill = unlist(paper_Transcriptome.Subtype))) +
  geom_bar(position = "fill",
           color = "black") +
  theme_classic() +
  xlab("Group") +
  ylab("Proportion") +
  scale_fill_manual(name = "Subtype", 
                    values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "snow2"),
                    labels = c("Classical","Mesenchymal", "Neural", "Proneural", "Unclassified")) + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  xlab("")

pdf("results/TCGA/Subtypes_across_PRNP-High_and_PRNP-Low_groups.pdf",
    width = 4,
    height = 4.5)
groups.prnp.sub
dev.off()

groups.prnp.exp.hist <- ggplot(PRNP_counts,
       aes(x = ENSG00000171867.17,
           y = ..density..,
           fill = PRNP_status)) +
  geom_histogram(color = "black",
                 alpha=0.6,
                 position = "identity",
                 binwidth = 0.025) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  ylab("Density") +
  xlab("PRNP expression") +
  scale_fill_manual(name = "Group",
                    values = c("orange", "steelblue1")) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))

pdf("results/TCGA/PRNP_distribution_PRNP-High_and_PRNP-Low_groups.pdf",
    width = 7,
    height = 3)
groups.prnp.exp.hist
dev.off()



genesets <- read.csv("data/c5.go.v7.5.1.symbols.gmt",
                     sep = "\t",
                     header = FALSE,
                     row.names = 1)
genesets <- as.data.frame(t(genesets))


gene_survival <- function(gene) {
  gene.exp <- pGBM_log_metadata[pGBM_log_metadata$paper_IDH.status == "Mutant",] %>% dplyr::select(gene, paper_Survival..months., paper_Vital.status..1.dead.)
  colnames(gene.exp) <- c("expression","months","status")
  cutoff.test <- maxstat.test(Surv(months, status) ~ expression, 
                              data = gene.exp,
                              smethod = "LogRank",
                              pmethod="HL")
  gene.exp <- gene.exp %>% mutate(classification = case_when(expression <= as.numeric(cutoff.test[["estimate"]]) ~ "Low",
                                                             expression > as.numeric(cutoff.test[["estimate"]]) ~ "High"))
  fit <- survfit(Surv(months, status) ~ classification, 
                 data = gene.exp,
                 type = "kaplan-meier")
  plot <- ggsurvplot(fit, 
                     data = gene.exp,
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = FALSE,
                     pval.coord = c(25, 0.75),
                     pval.method.coord = c(25, 0.85), 
                     ggtheme = theme_classic(),
                     legend = "right",
                     legend.labs = c(paste0(gene, " high"), paste0(gene, " low")),
                     xlab = "Time in months",
                     palette = c("tan1", "steelblue1"),
                     legend.title = paste0("Cutoff = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4)))
  plot <- ggpar(plot, 
                font.legend = 15,
                font.tickslab = 15,
                font.x = 15,
                font.y = 15)
  return(plot)
}


gene_survival("ENSG00000171867.17")


cowplot::plot_grid(plotlist = list(groups.prnp.exp.hist, 
                                   NULL,
                                   groups.prnp.idh, 
                                   groups.prnp.sub),
                   align = "hv",
                   ncol = 2,
                   nrow = 2,
                   axis = "btlr")

########## Step 9: Identifying differentially expressed transcripts ##########

dds <- DESeqDataSetFromMatrix(countData = PRNP_quartiles_count_data, 
                              colData = PRNP_quartiles_status,
                              design = ~ PRNP_status)

dds <- dds[rowSums(counts(dds)) > 1,]

# Estimando size factors.
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

res_unshrunken <- results(dds, 
                          contrast=c('PRNP_status','PRNP-High', 'PRNP-Low'), 
                          alpha = 0.05)

res_shrunken <- lfcShrink(dds = dds, 
                          res = res_unshrunken, 
                          contrast=c('PRNP_status','PRNP-High', 'PRNP-Low'), 
                          type = "normal")

deg.sig <- subset(res_shrunken, padj < 0.05)

print("Identified differentally expressed transcripts between PRNP-High and PRNP-Low GBM samples.")

deg.sig <- as.data.frame(deg.sig)

deg.all <- as.data.frame(res_shrunken)
deg.all <- deg.all %>% mutate(classification = case_when(padj <= 0.05 & log2FoldChange > 0 ~ "Upregulated",
                                                         padj <= 0.05 & log2FoldChange < 0 ~ "Downregulated",
                                                         padj > 0.05 ~ "Non-significant",
                                                         is.na(padj) ~ "Non-significant"))

ggplot(deg.all,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = classification)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_manual(values = c("blue", "lightgray", "red"),
                     name = "",
                     labels = c(paste0("Upregulated (n=", nrow(deg.all[deg.all$classification == "Upregulated",]), ")"),
                     paste0("Non-significant (n=", nrow(deg.all[deg.all$classification == "Non-significant",]), ")"),
                     paste0("Downregulated (n=", nrow(deg.all[deg.all$classification == "Downregulated",]), ")"))) +
  xlab(expression(log[2]("Fold change"))) +
  ylab(expression(-log[10]("Adjusted p-value"))) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.text.align = 0) + 
  guides(color = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

pdf("results/TCGA/volcano.pdf",
    width = 6,
    height = 4)
volcano
dev.off()

pdf("results/TCGA/arranged_plots_groups.pdf",
    width = 20,
    height = 3.5)
groups <- cowplot::plot_grid(plotlist = list(groups.prnp.exp, groups.prnp.exp.hist, 
                                   groups.prnp.idh, groups.prnp.sub, volcano),
                   align = "hv",
                   ncol = 5,
                   nrow = 1,
                   axis = "btlr",
                   rel_widths = c(2.5,4,1.6,1.6,2.5))
dev.off()

write.csv(deg.sig,
          file = "results/TCGA/TCGA-GBM_DETs_padj005.csv")

pdf("results/TCGA/arranged_plots.pdf",
    width = 20,
    height = 8)
cowplot::plot_grid(plotlist = list(all, groups),
                   align = "hv",
                   ncol = 1,
                   nrow = 2,
                   axis = "btrl")
dev.off()

#################################3

PRNP.median <- median(PRNP_counts_all$PRNP)

PRNP_counts_km <- PRNP_counts_all %>%
  mutate(PRNP_status_median = case_when(PRNP_counts_all$PRNP <= PRNP.median ~ "PRNP low",
                                      PRNP_counts_all$PRNP > PRNP.median ~ "PRNP high"))

PRNP_counts_km <- cbind(PRNP_counts_km, GBM_primary_tumor_metadata)

km_fit <- survfit(Surv(unlist(paper_Survival..months.), unlist(paper_Vital.status..1.dead.)) ~ PRNP_status_median, 
                  data = PRNP_counts_km)

ggsurvplot(
  km_fit,                     # survfit object with calculated statistics.
  data = PRNP_counts_km,
  legend.labs = c("PRNP-High", "PRNP-Low"),
  pval = TRUE
)










########## Step 10: Finding correlation between PRNP and traffic signatures ##########

genesets <- read.csv("data/c5.go.v7.5.1.symbols.gmt",
                     header = FALSE,
                     row.names = 1,
                     sep = "\t")
genesets.t <- t(genesets)
genesets.t <- as.data.frame(genesets.t)

traffic.genesets <- genesets.t %>% dplyr::select(grep(pattern = "VESICLE", 
                                                      x = colnames(genesets.t), 
                                                      value = TRUE),
                                                 GOBP_EXOCYTOSIS,
                                                 GOBP_ENDOCYTOSIS,
                                                 GOBP_SECRETION,
                                                 GOBP_INTRACELLULAR_TRANSPORT,
                                                 GOBP_INTRACELLULAR_PROTEIN_TRANSPORT,
                                                 GOBP_TRANSPORT_ALONG_MICROTUBULE)

genesets.list <- list()
for (i in colnames(traffic.genesets)) {
  genesets.list[[tolower(i)]] <- traffic.genesets[,i][traffic.genesets[,i] != ""]
}

traffic.gsva <- gsva(as.matrix(pGBM_log),
                     genesets.list)
rownames(traffic.gsva) <- gsub(x = rownames(traffic.gsva),
                               pattern = c("gobp_","gocc_","gomf_","hp_"),
                               replacement = "")

traffic.gsva <- as.data.frame(t(traffic.gsva))

expr <- cbind(traffic.gsva, t(pGBM_log)[,"PRNP"])
names(expr)[names(expr) == 't(pGBM_log)[, "PRNP"]'] <- "PRNP"

corr <- rcorr(as.matrix(expr),
              type = "spearman")
corr.prnp <- cbind(corr$r[,"PRNP"], corr$P[,"PRNP"])
colnames(corr.prnp) <- c("R", "P")
corr.prnp <- as.data.frame(corr.prnp)
corr.prnp <- corr.prnp[-c(126),]

ggplot(corr.prnp,
       aes(x = reorder(rownames(corr.prnp), R),
           fill = P,
           size = R)) +
  geom_dotplot() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 4)) +
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(3, "Reds")))(1000),
                       limits = c(0,0.05)) +
  coord_flip()




print("Finished analyzing GBM data from TCGA.")


library(gsoap)

data("pxgenes")

up <- read.csv("results/TCGA/TCGA_upDETs_gProfiler_ResultsList.csv",
               sep= ";",
               row.names = NULL,
               header = T)
down <- read.csv("results/TCGA/TCGA_downDETs_gProfiler_ResultsList.csv",
               sep= ";",
               row.names = NULL,
               header = T)

up$group <- "up"
down$group <- "down"

table <- rbind(up, down)

library(ggplot2)

ggplot(up[up$negative_log10_of_adjusted_p_value > 15,],
       aes(y = reorder(term_name, -log10(adjusted_p_value)),
           x = group,
           fill = -log10(adjusted_p_value))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


emapplot(up)

up$color <- ifelse(grepl("vesicle membrane|intracellular vesicle|transport vesicle|secretory vesicle|endocytic vesicle|vesicle-mediated transport|vesicle membrane|synaptic vesicle|extracellular exosome", up$term_name), "color", "no color")

tiff("results/TCGA/pathways.tif",
     units = "in",
     res = 600,
     width = 9,
     height = 4)
ggplot(up %>% arrange(color, .by_group = TRUE),
       aes(y = negative_log10_of_adjusted_p_value,
           x = intersection_size,
           label= term_name,
           color = color)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_manual(values = c("red","snow4")) +
  ylab(expression(-log[10]("Adjusted p-value"))) +
  xlab("Number of transcripts") +
  geom_text_repel(aes(label=ifelse(grepl("vesicle membrane|intracellular vesicle|transport vesicle|secretory vesicle|endocytic vesicle|vesicle-mediated transport|vesicle membrane|synaptic vesicle|extracellular exosome", 
                                         term_name),
                                   as.character(term_name),
                                   '')),
                  max.overlaps = Inf,
                  direction = "y",
                  xlim = c(900, NA),
                  ylim = c(2, NA),
                  nudge_x = 10,
                  nudge_y = 32, #23, 35
                  size = 4, 
                  force = 1,
                  hjust = 5,
                  seed = 1, 
                  segment.alpha = 0.5,
                  min.segment.length = 0,
                  segment.size = 1) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.position = "none") +
  ylim(c(0,80))
dev.off()    






