library(survminer)
library(ggplot2)
library(dplyr)
library(GSVA)
library(maxstat)
library(survival)

pGBM_log_metadata <- readRDS("results/TCGA/PrimaryGBMs_logData+Metadata.RDS")
genesets <- read.csv("data/c5.go.v7.5.1.symbols.gmt",
                     sep = "\t",
                     header = FALSE,
                     row.names = 1)
genesets <- as.data.frame(t(genesets))

pGBM_log_metadata <- pGBM_log_metadata[pGBM_log_metadata$paper_IDH.status == "WT",]

signature_survival <- function(geneset.name) {
  genes.list <- list(geneset.name = genesets[,geneset.name][genesets[,geneset.name] != ""])
  signature.scores <- gsva(t(pGBM_log_metadata[,106:56707]),
                 genes.list)
  signature.scores <- as.data.frame(t(signature.scores))
  signature.scores <- cbind(signature.scores, pGBM_log_metadata %>% select(paper_Survival..months., paper_Vital.status..1.dead.))
  colnames(signature.scores) <- c("scores","months","status")
  cutoff.test <- maxstat.test(Surv(months, status) ~ scores, 
                              data = signature.scores,
                              smethod = "LogRank",
                              pmethod = "HL")
  signature.scores <- signature.scores %>% mutate(classification = case_when(scores <= as.numeric(cutoff.test[["estimate"]]) ~ "Low",
                                                                             scores > as.numeric(cutoff.test[["estimate"]]) ~ "High"))
  fit <- survfit(Surv(months, status) ~ classification, 
                 data = signature.scores,
                 type = "kaplan-meier")
  dist <- ggplot(signature.scores,
         aes(x = "",
             y = scores)) +
    geom_boxplot() +
    geom_jitter(aes(color = classification)) +
    theme_classic() +
    ylab("GSVA scores") +
    xlab("") +
    theme(axis.ticks.x = element_blank(),
          legend.position = "top") +
    scale_color_discrete(name = "Signature")
  surv.plot <- ggsurvplot(fit, 
                     data = signature.scores,
                     pval = TRUE,
                     pval.method = TRUE,
                     conf.int = FALSE,
                     pval.coord = c(25, 0.75),
                     pval.method.coord = c(25, 0.85),
                     ggtheme = theme_classic(),
                     legend = "right",
                     legend.labs = c("Above cutoff", "Below cutoff"),
                     xlab = "Time in months",
                     title = geneset.name,
                     legend.title = paste0("Cutoff = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))) 
  plot.list <- list("Dist" = dist,
                    "Surv" = ggpubr::ggpar(surv.plot,
                                           font.legend = list(size = 15),
                                           font.x = list(size = 15),
                                           font.y = list(size = 15),
                                           font.xtickslab = list(size = 15),
                                           font.ytickslab = list(size = 15)))
  return(plot.list)
}


p1 <- signature_survival("GOCC_VESICLE_LUMEN") 
p2 <- signature_survival("GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS")
p3 <- signature_survival("GOBP_EXOCYTOSIS")
p4 <- signature_survival("GOBP_INTERCELLULAR_TRANSPORT")
p5 <- signature_survival("GOBP_EXPORT_FROM_CELL")
p6 <- signature_survival("GOCC_EXOCYTIC_VESICLE")
p7 <- signature_survival("GOBP_REGULATION_OF_VESICLE_FUSION")
p8 <- signature_survival("GOBP_EXOCYTIC_PROCESS")
p9 <- signature_survival("GOBP_VESICLE_DOCKING")
p10 <- signature_survival("GOCC_TRANSPORT_VESICLE")

pdf("surv_plots.pdf",
    width = 10,
    height = 10)
plots <- arrange_ggsurvplots(x = list(p1[[2]],
                                   p2[[2]],
                                   p3[[2]],
                                   p4[[2]],
                                   p5[[2]],
                                   p6[[2]],
                                   p7[[2]],
                                   p8[[2]],
                                   p9[[2]],
                                   p10[[2]]),
                   ncol = 2,
                   nrow = 5,
                   print = FALSE)
ggsave(filename = "surv_plots.pdf", plots,
       width = 12,
       height = 20)



pdf("results/TCGA/GOCC_VESICLE_LUMEN_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOCC_VESICLE_LUMEN_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p1[[2]]
dev.off()


p <- signature_survival("GOBP_REGULATION_OF_VESICLE_FUSION") 

pdf("results/TCGA/GOBP_REGULATION_OF_VESICLE_FUSION_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_REGULATION_OF_VESICLE_FUSION_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()

signature_survival("GOCC_SECRETORY_VESICLE") 
signature_survival("GOCC_VESICLE_MEMBRANE") 
signature_survival("GOBP_ENDOCYTOSIS") 
signature_survival("GOBP_SECRETION")
signature_survival("GOBP_MEMBRANE_INVAGINATION") 
signature_survival("GOBP_VESICLE_ORGANIZATION")
signature_survival("GOCC_SECRETORY_VESICLE")

p <- signature_survival("GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS")

pdf("results/TCGA/GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()

p <- signature_survival("GOBP_EXOCYTOSIS")

pdf("results/TCGA/GOBP_EXOCYTOSIS_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_EXOCYTOSIS_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()

p <- signature_survival("GOBP_EXOCYTIC_PROCESS")

pdf("results/TCGA/GOBP_EXOCYTIC_PROCESS_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_EXOCYTIC_PROCESS_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()

p <- signature_survival("GOBP_INTERCELLULAR_TRANSPORT")

pdf("results/TCGA/GOBP_INTERCELLULAR_TRANSPORT_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_INTERCELLULAR_TRANSPORT_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()

p <- signature_survival("GOBP_EXPORT_FROM_CELL")

pdf("results/TCGA/GOBP_EXPORT_FROM_CELL_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOBP_EXPORT_FROM_CELL_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()


p <- signature_survival("GOCC_EXOCYTIC_VESICLE")

pdf("results/TCGA/GOCC_EXOCYTIC_VESICLE_distribution_optimal_cutoff.pdf",
    width = 3,
    height = 4)
p[[1]]
dev.off()

pdf("results/TCGA/GOCC_EXOCYTIC_VESICLE_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
p[[2]]
dev.off()



library(ggpubr)
arrange_ggsurvplots(x = list(p1,p2,p3,p4,p5,p6,p7,p8),
          align = "hv",
          ncol = 4,
          nrow = 2)

p1+p2


gene_survival <- function(gene) {
  gene.exp <- pGBM_log_metadata %>% select(gene, paper_Survival..months., paper_Vital.status..1.dead.)
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

pdf("results/TCGA/PRNP_survival_optimal_cutoff.pdf",
    width = 6,
    height = 4)
gene_survival("PRNP")
dev.off()

pdf("results/TCGA/DSTN_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
gene_survival("DSTN")
dev.off()

pdf("results/TCGA/ANXA1_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
gene_survival("ANXA1")
dev.off()

pdf("results/TCGA/RAB31_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
gene_survival("RAB31")
dev.off()

pdf("results/TCGA/SYPL1_survival_optimal_cutoff.pdf",
    width = 5,
    height = 4)
gene_survival("SYPL1")
dev.off()