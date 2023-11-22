library(Hmisc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

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

tissue.log.metadata <- readRDS("results/TCGA-GBM/Tissue_logData+Metadata.RDS")

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

pdf("results/TCGA/42_genes_across_tissues.pdf",
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


pGBM.log.metadata <- readRDS("results/TCGA/PrimaryGBMs_logData+Metadata.RDS")

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


pdf("results/TCGA/42_genes_across_idh.pdf",
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

pdf("results/TCGA/42_genes_across_subtype.pdf",
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

