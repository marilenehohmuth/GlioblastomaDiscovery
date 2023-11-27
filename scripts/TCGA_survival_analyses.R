
# Load packages -----------------------------------------------------------

library(survminer)                # 0.4.9
library(ggplot2)                  # 3.4.4
library(dplyr)                    # 1.1.3
library(GSVA)                     # 1.48.2
library(maxstat)                  # 0.7-25
library(survival)                 # 3.5-7
library(org.Hs.eg.db)             # 3.17.0



# Load data ---------------------------------------------------------------

# Load primary GBM data.
pGBM_log_metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/survival_analysis/PrimaryGBMs_logData+Metadata_filtered.RDS"))
pGBM_log <- pGBM_log_metadata[,grepl("ENSG", colnames(pGBM_log_metadata))] # Only genes!
colnames(pGBM_log) <- gsub("\\..*", "", colnames(pGBM_log))

# Load GO gene sets.
genesets <- read.csv(
    paste0(getwd(), "/data/c5.go.v7.5.1.symbols.gmt"),
    sep = "\t",
    header = FALSE,
    row.names = 1)

# Remove column that contains accession links & transpose dataframe.
genesets$V2 <- NULL 
genesets <- as.data.frame(t(genesets))

# Signature survival ------------------------------------------------------

signature_survival <- function(geneset.name, gene) {

    # Get genes associated with user-provided gene set name.
    genes.list <- list(geneset.name = genesets[,geneset.name][genesets[,geneset.name] != ""])
    
    # Get ENSEMBL IDs.
    correspondence <- AnnotationDbi::select(
        org.Hs.eg.db, 
        keys = unlist(genes.list), 
        columns = "ENSEMBL", 
        keytype = "SYMBOL"
    )
    genes.list <- list(correspondence$ENSEMBL)

    # Calculate Gene Set Variation Analysis (GSVA) scores.
    signature.scores <- gsva(t(pGBM_log), genes.list)
    signature.scores <- as.data.frame(t(signature.scores))

    # Concatenate GSVA scores and survival information.
    signature.scores <- cbind(signature.scores, pGBM_log_metadata %>% dplyr::select(paper_Survival..months., paper_Vital.status..1.dead.))
    colnames(signature.scores) <- c("scores", "months", "status")
    
    # Get optimal cut-off.
    cutoff.test <- maxstat.test(
        Surv(months, status) ~ scores, 
        data = signature.scores,
        smethod = "LogRank",
        pmethod = "HL"
    )

    # Add sample classification based on given cut-off.
    signature.scores <- signature.scores %>% mutate(classification = case_when(
        scores <= as.numeric(cutoff.test[["estimate"]]) ~ "Low",
        scores > as.numeric(cutoff.test[["estimate"]]) ~ "High"
    ))

    # Carry out survival estimate.
    fit <- survfit(
        Surv(months, status) ~ classification, 
        data = signature.scores,
        type = "kaplan-meier"
    )

    # Create survival plot with risk table.
    surv.plot <- ggsurvplot(
        fit, 
        data = signature.scores,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = FALSE,
        pval.coord = c(25, 0.75),
        pval.method.coord = c(25, 0.85),
        ggtheme = theme_classic(),
        legend = "right",
        legend.labs = c(paste0("Above cut-off (n=", nrow(signature.scores[signature.scores$classification == "High",]), ")"), paste0("Below cut-off (n=", nrow(signature.scores[signature.scores$classification == "Low",]), ")")),
        xlab = "Time in months",
        title = geneset.name,
        risk.table = TRUE,
        legend.title = paste0("Cut-off = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))
    )

    # surv.plot[[1]] = survival curves.
    # surv.plot[[2]] = risk table.
    # surv.plot[[3]] = dataframe used for survival analysis.

    # Adjust text size of survival plot.
    i <- 1
    while (i < 3) {
        surv.plot[[i]] <- surv.plot[[i]] + 
            theme(
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 15),
                legend.title = element_text(size =15)
            )
        i <- i + 1
    }

    # Save survival plot with risk table.
    pdf(paste0(getwd(), "/results/TCGA-GBM/survival_analysis/", geneset.name, "_survival_analysis.pdf"), width = 10, height = 6)
    print(surv.plot)
    dev.off()

    # Create plot showing the expression of a given gene in the high and low groups.
    correspondence <- AnnotationDbi::select(
        org.Hs.eg.db, 
        keys = gene, 
        columns = "ENSEMBL", 
        keytype = "SYMBOL"
    )
    gene.exp <- pGBM_log %>% dplyr::select(correspondence$ENSEMBL)
    signature.scores <- cbind(signature.scores, gene.exp)
    colnames(signature.scores) <- c("scores", "months", "status", "classification", "gene")
    gene.plot <- ggplot(
        signature.scores,
        aes(x = classification, y = as.numeric(gene))
    ) +
        geom_violin(scale = "width", aes(fill = classification)) +
        geom_boxplot(width = 0.25, outlier.shape = NA) +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12),
            legend.position = "none",
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        ) +
        ylab(paste0(gene, " expression")) +
        scale_x_discrete(
            labels = c(paste0("Above cut-off (n=", nrow(signature.scores[signature.scores$classification == "High",]), ")"), paste0("Below cut-off (n=", nrow(signature.scores[signature.scores$classification == "Low",]), ")"))
        ) +
        ggtitle(geneset.name) +
        stat_compare_means(
            label = "p.format",
            label.y.npc = "bottom",
            label.x.npc = "left",
            size = 4,
            method = "wilcox.test"
        ) 
    
    pdf(paste0(getwd(), "/results/TCGA-GBM/survival_analysis/", geneset.name, "_", gene, "_expression_in_groups.pdf"), width = 3, height = 3.5)
    print(gene.plot)
    dev.off()
}


## Calculate signature -----------------------------------------------------

signature_survival("GOCC_VESICLE_LUMEN", "PRNP") 
signature_survival("GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS", "PRNP")
signature_survival("GOBP_EXOCYTOSIS", "PRNP")
signature_survival("GOBP_INTERCELLULAR_TRANSPORT", "PRNP")
signature_survival("GOBP_EXPORT_FROM_CELL", "PRNP")
signature_survival("GOCC_EXOCYTIC_VESICLE", "PRNP")
signature_survival("GOBP_REGULATION_OF_VESICLE_FUSION", "PRNP")
signature_survival("GOBP_EXOCYTIC_PROCESS", "PRNP")
signature_survival("GOBP_VESICLE_DOCKING", "PRNP")
signature_survival("GOCC_TRANSPORT_VESICLE", "PRNP")
signature_survival("GOCC_SECRETORY_VESICLE", "PRNP") 
signature_survival("GOCC_VESICLE_MEMBRANE", "PRNP") 
signature_survival("GOBP_ENDOCYTOSIS", "PRNP") 
signature_survival("GOBP_SECRETION", "PRNP")
signature_survival("GOBP_MEMBRANE_INVAGINATION", "PRNP") 
signature_survival("GOBP_VESICLE_ORGANIZATION", "PRNP")

# Gene survival ------------------------------------------------------

gene_survival <- function(gene) {

    correspondence <- AnnotationDbi::select(
        org.Hs.eg.db, 
        keys = gene, 
        columns = "ENSEMBL", 
        keytype = "SYMBOL"
    )

    # Get gene expression and survival data.
    gene.exp <- pGBM_log %>% dplyr::select(correspondence$ENSEMBL)
    gene.exp <- cbind(gene.exp, pGBM_log_metadata %>% dplyr::select(paper_Survival..months., paper_Vital.status..1.dead.))
    colnames(gene.exp) <- c("expression", "months", "status")
    
    # Calcualate optimal cut-off.
    cutoff.test <- maxstat.test(
        Surv(months, status) ~ expression, 
        data = gene.exp,
        smethod = "LogRank",
        pmethod="HL"
    )
  
    # Add classification to samples based on given cut-off.
    gene.exp <- gene.exp %>% mutate(classification = case_when(
        expression <= as.numeric(cutoff.test[["estimate"]]) ~ "Low",
        expression > as.numeric(cutoff.test[["estimate"]]) ~ "High")
    )

    # Carry out survival estimate.
    fit <- survfit(Surv(months, status) ~ classification, 
                    data = gene.exp,
                    type = "kaplan-meier")
  
    surv.plot <- ggsurvplot(
        fit, 
        data = gene.exp,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = FALSE,
        pval.coord = c(25, 0.75),
        pval.method.coord = c(25, 0.85),
        ggtheme = theme_classic(),
        legend = "right",
        legend.labs = c(paste0("Above cut-off (n=", nrow(gene.exp[gene.exp$classification == "High",]), ")"), paste0("Below cut-off (n=", nrow(gene.exp[gene.exp$classification == "Low",]), ")")),
        xlab = "Time in months",
        title = gene,
        risk.table = TRUE,
        legend.title = paste0("Cut-off = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))
    )

    # surv.plot[[1]] = survival curves.
    # surv.plot[[2]] = risk table.
    # surv.plot[[3]] = dataframe used for survival analysis.

    # Adjust text size of survival plot.
    i <- 1
    while (i < 3) {
        surv.plot[[i]] <- surv.plot[[i]] + 
            theme(
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 15),
                legend.title = element_text(size =15),
                plot.title = element_text(face = "italic")
            )
        i <- i + 1
    }

    pdf(paste0(getwd(), "/results/TCGA-GBM/survival_analysis/", gene, "_survival_analysis.pdf"), width = 10, height = 6)
    print(surv.plot)
    dev.off()
}

genes <- read.csv(
    paste0(getwd(), "/results/comparison_all/all_datasets_common_PRNPhigh_or_positive_upGenes.csv"),
    header = TRUE,
    row.names = 1
)

for(gene in genes$x) {
    gene_survival(gene)
}