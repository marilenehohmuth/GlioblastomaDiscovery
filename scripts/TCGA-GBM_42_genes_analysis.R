<<<<<<< HEAD:scripts/TCGA-GBM_43_genes_analysis.R
library(Hmisc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# Load PRNP-High vs PRNP-Low DEGs (TCGA-GBM, bulk RNA-seq)
tcga <- read.csv(
  paste0(getwd(), "/results/TCGA-GBM/TCGA-GBM_PRNP-High_vs_PRNP-Low_DETs_padj0.05.csv"),
  header = TRUE,
  row.names = 1
)
# Get upregulated genes in the PRNP-High samples.
tcga <- tcga$gene_symbol[tcga$classification == "Upregulated"] 

# Load PRNP+ and PRNP- marker genes (Darmanis et al, single-cell RNA-seq)
darmanis <- read.csv(
  paste0(getwd(), "/results/darmanis/Darmanis_PRNP+_vs_PRNP-_signature_marker_genes.csv"),
  header = TRUE,
  row.names = 1
)
# Get upregulated genes in the PRNP+ cell population.
darmanis <- darmanis$gene[darmanis$population == "PRNP_positive_cells" & darmanis$fc > 0 & darmanis$pval <= 0.05]

# Load PRNP+ and PRNP- marker genes (Neftel et al, single-cell RNA-seq)
neftel <- read.csv(
  "/nfs/team205/jb62/other/GlioblastomaDiscovery/results/neftel/Neftel_PRNP+_vs_PRNP-_signature_marker_genes.csv",
  header = TRUE,
  row.names = 1
)
# Get upregulated genes in the PRNP+ cell population.
neftel <- neftel$gene[neftel$population == "PRNP_positive_cells" & neftel$fc > 0 & neftel$pval <= 0.05]

# Load PRNP+ and PRNP- marker genes (Richards et al, single-cell RNA-seq)
richards <- read.csv(
  "/nfs/team205/jb62/other/GlioblastomaDiscovery/results/richards/Richards_PRNP+_vs_PRNP-_signature_marker_genes.csv",
  header = TRUE,
  row.names = 1
)
# Get upregulated genes in the PRNP+ cell population.
richards <- richards$gene[richards$population == "PRNP_positive_cells" & richards$fc > 0 & richards$pval <= 0.05]

# Create list with gene lists.
gene_lists <- list(
  "TCGA" = tcga,
  "Darmanis" = darmanis,
  "Neftel" = neftel,
  "Richards" = richards
)

# Create Venn diagram.
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

# Get genes common to all 4 datasets.
common_upGenes <- Reduce(intersect, gene_lists)

# Get ENSEMBL IDs.
correspondence <- select(
  org.Hs.eg.db, 
  keys = common_upGenes, 
  columns = "ENSEMBL", 
  keytype = "SYMBOL"
)

# Load TCGA-GBM data.
tissue.log.metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/Tissue_logData+Metadata.RDS"))

# Removing version identifier from ENSEMBL IDs.
colnames(tissue.log.metadata) <- gsub("\\..*", "", colnames(tissue.log.metadata))

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

df_tissues$symbol <- correspondence$SYMBOL[match(df_tissues$gene, correspondence$ENSEMBL)]

pdf(
  paste0(getwd(), "/results/comparison_all/all_datasets_common_genes_across_tissues.pdf"),
  width = 14,
  height = 20
)
ggplot(
  df_tissues[df_tissues$symbol != "PRNP",],
  aes(x = tissue,  y = as.numeric(expression))
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
    strip.text.x = element_text(face = "italic")) +
  stat_compare_means( # Kruskal-Wallis test.
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "center",
    size = 3,
    method = 'kruskal.test'
  )
dev.off()


pGBM.log.metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/PrimaryGBMs_logData+Metadata.RDS"))

# Removing version identifier from ENSEMBL IDs.
colnames(pGBM.log.metadata) <- gsub("\\..*", "", colnames(pGBM.log.metadata))

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

df_idh$symbol <- correspondence$SYMBOL[match(df_idh$gene, correspondence$ENSEMBL)]

pdf(
  paste0(getwd(), "/results/comparison_all/all_datasets_common_genes_across_idh.pdf"),
  width = 14,
  height = 20
)
ggplot(
  df_idh,
  aes(x = idh, y = as.numeric(expression))) +
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

df_subtype$symbol <- correspondence$SYMBOL[match(df_subtype$gene, correspondence$ENSEMBL)]

pdf(
  paste0(getwd(), "/results/comparison_all/all_datasets_common_genes_across_subtypes.pdf"),
  width = 14,
  height = 20
)
ggplot(
  df_subtype[df_subtype$symbol != "PRNP",],
  aes(x = subtype, y = as.numeric(expression))) +
  geom_violin(aes(fill = subtype), scale = "width") +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
  facet_wrap(~symbol, scales = "free") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "none",
    strip.text.x = element_text(face = "italic")) +
  stat_compare_means(
    label = "p.format",
    label.y.npc = "bottom",
    label.x.npc = "center",
    size = 3,
    method = 'kruskal.test'
  )+
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#7CAE00", "gray"))
dev.off()

=======
library(Hmisc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(biomaRt)                  # 2.46.3 

common.genes <- c('DSTN',
                  'PIGT',
                  'HADHB',
                  'ATP6AP2',
                  'FOSB',
                  'PEBP1',
                  'PLTP',
                  'ZFP36L2',
                  'DUSP1',
                  'ITM2C',
                  'YWHAB',
                  'CD59',
                  'CLU',
                  'C21orf62',
                  'GFAP',
                  'MATN2',
                  'PON2',
                  'SPARC',
                  'F3',
                  'GATM',
                  'NDFIP1',
                  'CBR1',
                  'ITM2B',
                  'CST3',
                  'ZFP36',
                  'TIMP1',
                  'PTTG1IP',
                  'SBDS',
                  'SYPL1',
                  'TPST1',
                  'PMP22',
                  'CALM1',
                  'WLS',
                  'ANXA1',
                  'GJA1',
                  'SMOX',
                  'RHBDD2',
                  'HIGD1A',
                  'CLDND1',
                  'RAB31',
                  'REEP5',
                  'SERPINB6')

## load metadata
tissue.log.metadata <- readRDS("results/TCGA-GBM/Tissue_logData+Metadata.RDS")
## get ids and remove . after it
gene_ids <- names(tissue.log.metadata)[109:length(names(tissue.log.metadata))]
gene_ids <- sub("\\..*", "", gene_ids)
names(tissue.log.metadata)[109:length(names(tissue.log.metadata))] <- gene_ids

## change gene ids to symbols
## determine biomaRt settings.
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",  mart = ensembl)

## get correspondence between ENSEMBL IDs and gene symbols.
genemap <- getBM(
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = gene_ids, 
  mart = ensembl
)

## rename columns
for (i in seq_len(nrow(genemap))) {
  names(tissue.log.metadata)[names(tissue.log.metadata) == genemap$ensembl_gene_id[i]] <- genemap$hgnc_symbol[i]
}

df_tissues <- data.frame(matrix(ncol = 3))
for (gene in common.genes) {
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

pdf("results/TCGA-GBM/42_genes_across_tissues.pdf",
    width = 14,
    height = 16)
ggplot(df_tissues,
       aes(x = tissue,
           y = as.numeric(expression),
           fill = tissue)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        axis.text.y = element_text(size = 12)) +
  stat_compare_means(label = "p.format",
                     label.y.npc = "bottom",
                     label.x.npc = "center",
                     size = 4)+
  ylab(NULL) +
  xlab(NULL)
dev.off()


pGBM.log.metadata <- readRDS("results/TCGA-GBM/PrimaryGBMs_logData+Metadata.RDS")

## get ids and remove . after it
gene_ids <- names(pGBM.log.metadata)[110:length(names(pGBM.log.metadata))]
gene_ids <- sub("\\..*", "", gene_ids)
names(pGBM.log.metadata)[110:length(names(pGBM.log.metadata))] <- gene_ids

## change gene ids to symbols
## determine biomaRt settings.
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",  mart = ensembl)

## get correspondence between ENSEMBL IDs and gene symbols.
genemap <- getBM(
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = gene_ids, 
  mart = ensembl
)

## rename columns
for (i in seq_len(nrow(genemap))) {
  names(pGBM.log.metadata)[names(pGBM.log.metadata) == genemap$ensembl_gene_id[i]] <- genemap$hgnc_symbol[i]
}


df_idh <- data.frame(matrix(ncol = 3))
for (gene in common.genes) {
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_IDH.status == "WT"]) {
    df_idh <- rbind(df_idh, c(gene, expression, "IDH-WT"))
  }
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_IDH.status == "Mutant"]) {
    df_idh <- rbind(df_idh, c(gene, expression, "IDH-mutant"))
  }
  for (expression in pGBM.log.metadata[,gene][is.na(pGBM.log.metadata$paper_IDH.status) == TRUE]) {
    df_idh <- rbind(df_idh, c(gene, expression, "Unclassified"))
  }

}
df_idh <- df_idh[complete.cases(df_idh),]
colnames(df_idh) <- c("gene", "expression", "idh")
df_idh$idh <- factor(df_idh$idh, c("IDH-mutant", "IDH-WT", "Unclassified"))


pdf("results/TCGA-GBM/42_genes_across_idh.pdf",
    width = 14,
    height = 16)
ggplot(df_idh,
       aes(x = idh,
           y = as.numeric(expression),
           fill = idh)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        axis.text.y = element_text(size = 12)) +
  stat_compare_means(label = "p.format",
                     label.y.npc = "bottom",
                     label.x.npc = "center",
                     size = 4)+
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "gray"))
dev.off()



df_subtype <- data.frame(matrix(ncol = 3))
for (gene in common.genes) {
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome.Subtype == "CL"]) {
    df_subtype <- rbind(df_subtype, c(gene, expression, "Classical"))
  }
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome.Subtype == "NE"]) {
    df_subtype <- rbind(df_subtype, c(gene, expression, "Neural"))
  }
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome.Subtype == "PN"]) {
    df_subtype <- rbind(df_subtype, c(gene, expression, "Proneural"))
  }
  for (expression in pGBM.log.metadata[,gene][pGBM.log.metadata$paper_Transcriptome.Subtype == "ME"]) {
    df_subtype <- rbind(df_subtype, c(gene, expression, "Mesenchymal"))
  }
  for (expression in pGBM.log.metadata[,gene][is.na(pGBM.log.metadata$paper_Transcriptome.Subtype) == TRUE]) {
    df_subtype <- rbind(df_subtype, c(gene, expression, "Unclassified"))
  }
  
}
df_subtype <- df_subtype[complete.cases(df_subtype),]
colnames(df_subtype) <- c("gene", "expression", "subtype")
df_subtype$subtype <- factor(df_subtype$subtype, c("Proneural", "Classical", "Mesenchymal", "Neural", "Unclassified"))

pdf("results/TCGA-GBM/42_genes_across_subtype.pdf",
    width = 14,
    height = 16)
ggplot(df_subtype,
       aes(x = subtype,
           y = as.numeric(expression),
           fill = subtype)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25,
               outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        axis.text.y = element_text(size = 12)) +
  stat_compare_means(label = "p.format",
                     label.y.npc = "bottom",
                     label.x.npc = "center",
                     size = 4)+
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#7CAE00", "#00BFC4", "gray"))
dev.off()

>>>>>>> 6193e083941f03d8eb9303e56888490ce249cce4:scripts/TCGA-GBM_42_genes_analysis.R
