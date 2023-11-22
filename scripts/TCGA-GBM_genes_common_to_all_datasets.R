# TCGA-GBM ANALYSIS #########################################################
# @ This script generates the plots present in the third main figure of the # 
# @ manuscript.                                                             #
#############################################################################


##########################
#### Loading packages ####
##########################

# R 4.3.2
library(Hmisc)                    # 5.1_1
library(ggplot2)                  # 3.4.4
library(RColorBrewer)             # 1.1_3
library(ggpubr)                   # 0.6.0
library(VennDiagram)              # 1.7.3
library(dplyr)                    # 1.1.3
library(clusterProfiler)          # 4.8.1
library(org.Hs.eg.db)             # 3.17.0

################################################
#### Loading results from all four datasets ####
################################################

# Load PRNP-High vs PRNP-Low DEGs (TCGA-GBM, bulk RNA-seq)
tcga <- read.csv(
    paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_DETs_padj0.05.csv"),
    header = TRUE,
    row.names = 1
)
# Get upregulated genes in the PRNP-High samples.
tcga <- tcga$gene_symbol[tcga$classification == "Upregulated"]


# Load PRNP+ and PRNP- marker genes (Darmanis et al, single-cell RNA-seq)
darmanis <- read.csv(
    paste0(getwd(), "/results/darmanis/Darmanis_PRNP+_vs_PRNP-_signature_marker_genes.csv"),
    header = TRUE,
    row.names = NULL
)
# Get upregulated genes in the PRNP+ cell population.
darmanis <- darmanis$gene[darmanis$population == "PRNP_positive_cells" & darmanis$fc > 0 & darmanis$pval <= 0.05]


# Load PRNP+ and PRNP- marker genes (Neftel et al, single-cell RNA-seq)
neftel <- read.csv(
    "/nfs/team205/jb62/other/GlioblastomaDiscovery/results/neftel/Neftel_PRNP+_vs_PRNP-_signature_marker_genes.csv",
    header = TRUE,
    row.names = NULL
)
# Get upregulated genes in the PRNP+ cell population.
neftel <- neftel$gene[neftel$population == "PRNP_positive_cells" & neftel$fc > 0 & neftel$pval <= 0.05]


# Load PRNP+ and PRNP- marker genes (Richards et al, single-cell RNA-seq)
richards <- read.csv(
    "/nfs/team205/jb62/other/GlioblastomaDiscovery/results/richards/Richards_PRNP+_vs_PRNP-_signature_marker_genes.csv",
    header = TRUE,
    row.names = NULL
)
# Get upregulated genes in the PRNP+ cell population.
richards <- richards$gene[richards$population == "PRNP_positive_cells" & richards$fc > 0 & richards$pval <= 0.05]


# Create list with gene lists.
gene_lists <- list(
    "TCGA" = tcga,
    "Darmanis" = darmanis,
    "Neftel" = neftel,
    "Richards" = richards
)


##################################################################
#### Investigating genes that are common to all four datasets ####
##################################################################

# Create Venn diagram.
venn.diagram(
    x = gene_lists,
    category.names = c("TCGA", "Darmanis" , "Neftel" , "Richards"),
    filename = paste0(getwd(), '/results/comparison_all/all_datasets_common_genes_vennDiagram.png'),
    output = TRUE,
    imagetype = "png" ,
    height = 480 ,
    width = 480 ,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = brewer.pal(4, "Pastel1"),
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135, -135),
    cat.dist = c(0.055, 0.055, 0.085, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
)

# Get genes common to all 4 datasets.
common_upGenes <- Reduce(intersect, gene_lists)

# Get ENSEMBL IDs.
correspondence <- select(
    org.Hs.eg.db,
    keys = common_upGenes,
    columns = "ENSEMBL",
    keytype = "SYMBOL"
)

# Load TCGA-GBM data (including primary and recurrrent tumors & non-neoplastic tissue).
tissue.log.metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/Tissue_logData+Metadata.RDS"))

# Removing version identifier from ENSEMBL IDs.
colnames(tissue.log.metadata) <- gsub("\\..*", "", colnames(tissue.log.metadata))

# Create dataframe contaning each gene along with their expression in samples of each tissue type.
df_tissues <- data.frame(matrix(nrow = 0, ncol = 0))
for (gene in correspondence$ENSEMBL) {
    for (expression in tissue.log.metadata[,gene][tissue.log.metadata$definition == "Primary solid Tumor"]) {
        df_tissues <- rbind(df_tissues, c(gene, expression, "Primary"))
    }
    for (expression in tissue.log.metadata[,gene][tissue.log.metadata$definition == "Recurrent Solid Tumor"]) {
        df_tissues <- rbind(df_tissues, c(gene, expression, "Recurrent"))
    }
    for (expression in tissue.log.metadata[,gene][tissue.log.metadata$definition == "Solid Tissue Normal"]) {
        df_tissues <- rbind(df_tissues, c(gene, expression, "Non-neoplastic"))
    }
}
df_tissues <- df_tissues[complete.cases(df_tissues),]
colnames(df_tissues) <- c("gene", "expression", "tissue")

# Get gene symbols.
df_tissues$symbol <- correspondence$SYMBOL[match(df_tissues$gene, correspondence$ENSEMBL)]

# Generate plot & save it to output file.
pdf(
    paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_across_tissues.pdf"),
    width = 14,
    height = 20
)
ggplot(
    df_tissues[df_tissues$symbol != "PRNP",],
    aes(x = tissue, y = as.numeric(expression))
) +
    geom_violin(aes(fill = tissue), scale = "width") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    facet_wrap(~symbol, scales = "free") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(face = "italic")
    ) +
    stat_compare_means( # Kruskal-Wallis test.
        label = "p.format",
        label.y.npc = "bottom",
        label.x.npc = "center",
        size = 3,
        method = 'kruskal.test'
    )
dev.off()

# Load TCGA-GBM data (including only primary tumors).
pGBM.log.metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/PrimaryGBMs_logData+Metadata.RDS"))

# Removing version identifier from ENSEMBL IDs.
colnames(pGBM.log.metadata) <- gsub("\\..*", "", colnames(pGBM.log.metadata))

# Create dataframe contaning each gene along with their expression in samples of each IDH subtype.
df_idh <- data.frame(matrix(nrow = 0, ncol = 0))
for (gene in correspondence$ENSEMBL) {
    for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_IDH == "WT"]) {
        df_idh <- rbind(df_idh, c(gene, expression, "IDHwt"))
    }
    for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_IDH == "Mutant"]) {
    df_idh <- rbind(df_idh, c(gene, expression, "IDH-mutant"))
    }
    for (expression in pGBM.log.metadata[,gene][is.na(pGBM.log.metadata$paper_IDH)]) {
    df_idh <- rbind(df_idh, c(gene, expression, "Unclassified"))
    }
}
df_idh <- df_idh[complete.cases(df_idh),]
colnames(df_idh) <- c("gene", "expression", "idh")
df_idh$idh <- factor(df_idh$idh, c("IDH-mutant", "IDHwt", "Unclassified"))

# Get gene symbols.
df_idh$symbol <- correspondence$SYMBOL[match(df_idh$gene, correspondence$ENSEMBL)]

# Generate plot & save it to output file.
pdf(
    paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_across_idh.pdf"),
    width = 14,
    height = 20
)
ggplot(
    df_idh,
    aes(x = idh, y = as.numeric(expression))
) +
    geom_violin(aes(fill = idh), scale = "width") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    facet_wrap(~symbol, scales = "free") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(face = "italic")
    ) +
    stat_compare_means(
        label = "p.format",
        label.y.npc = "bottom",
        label.x.npc = "center",
        size = 3,
        method = 'kruskal.test'
    ) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D", "gray"))
dev.off()

# Create dataframe contaning each gene along with their expression in samples of each transcriptional subtype.
df_subtype <- data.frame(matrix(nrow = 0, ncol = 0))
for (gene in correspondence$ENSEMBL) {
    for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome == "CL"]) {
        df_subtype <- rbind(df_subtype, c(gene, expression, "Classical"))
    }
    for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome == "PN"]) {
        df_subtype <- rbind(df_subtype, c(gene, expression, "Proneural"))
    }
    for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome == "ME"]) {
        df_subtype <- rbind(df_subtype, c(gene, expression, "Mesenchymal"))
    }
    for (expression in pGBM.log.metadata[,gene][is.na(pGBM.log.metadata$paper_Transcriptome)]) {
        df_subtype <- rbind(df_subtype, c(gene, expression, "Unclassified"))
    }
}
df_subtype <- df_subtype[complete.cases(df_subtype),]
colnames(df_subtype) <- c("gene", "expression", "subtype")
df_subtype$subtype <- factor(df_subtype$subtype, c("Proneural", "Classical", "Mesenchymal", "Unclassified"))

# Get gene symbols.
df_subtype$symbol <- correspondence$SYMBOL[match(df_subtype$gene, correspondence$ENSEMBL)]

# Generate plot & save it to output file.
pdf(
    paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_across_subtypes.pdf"),
    width = 14,
    height = 20
)
ggplot(
    df_subtype[df_subtype$symbol != "PRNP",],
    aes(x = subtype, y = as.numeric(expression))
) +
    geom_violin(aes(fill = subtype), scale = "width") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    facet_wrap(~symbol, scales = "free") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(face = "italic")
    ) +
    stat_compare_means(
        label = "p.format",
        label.y.npc = "bottom",
        label.x.npc = "center",
        size = 3,
        method = 'kruskal.test'
    )+
    scale_fill_manual(values = c("#C77CFF", "#F8766D", "#7CAE00", "gray"))
dev.off()

# Perform Over-Representation Analysis (ORA).
comparison <- enrichGO(
    correspondence$SYMBOL,
    ont = "ALL",
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
)
comparison <- as.data.frame(comparison)

# Save ORA results to output file.
write.csv(
    comparison, 
    file = paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_ORA_enrichGO_ontALL_padj0.05.csv")
)

# Perform Over-Representation Analysis (ORA).
comparison_2 <- gost(
    query = correspondence$SYMBOL, 
    organism = "hsapiens", 
    ordered_query = FALSE, 
    multi_query = FALSE, 
    significant = TRUE, 
    exclude_iea = TRUE, 
    measure_underrepresentation = FALSE,
    evcodes = FALSE, 
    user_threshold = 0.05, 
    correction_method = "g_SCS", 
    domain_scope = "annotated", 
    custom_bg = NULL, 
    numeric_ns = "", 
    sources = NULL, 
    as_short_link = FALSE, 
    highlight = TRUE
)
comparison_2$result <- apply(comparison_2$result, 2, as.character)

# Save ORA results to output file.
write.table(
    as.data.frame(comparison_2$result), 
    file = paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_ORA_gProfiler2_ontALL_padj0.05.csv")
)

# Plot ORA results.
pdf(
    paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes_ORA_gProfiler2_padj0.05_allEnrichedTerms.pdf"),
    width = 8,
    height = 8
)
ggplot(
    as.data.frame(comparison_2$result),
    aes(x = -log10(as.numeric(p_value)), y = reorder(term_name, -log10(as.numeric(p_value))), size = as.numeric(intersection_size))
) +
    geom_point(shape = 21, color = "black", fill = "palegreen2") +
    theme_bw() +
    theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    xlab(expression(-log[10]~"(Adjusted p-value)")) +
    ylab("Term") +
    ggtitle("") +
    scale_size_continuous(range = c(1,5), name = "# Genes") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed")
dev.off()