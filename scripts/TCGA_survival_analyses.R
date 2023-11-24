
# Load packages -----------------------------------------------------------

library(survminer)                # 0.4-9
library(ggplot2)
library(dplyr)
library(GSVA)                     # 1.38.1
library(maxstat)                  # 0.7-25
library(survival)                 # 3.5-7
library(org.Hs.eg.db)             # 3.17.0



# Load data ---------------------------------------------------------------

# Load primary GBM data.
pGBM_log_metadata <- readRDS(paste0(getwd(), "/results/TCGA-GBM/PrimaryGBMs_logData+Metadata_filtered.RDS"))
pGBM_log <- pGBM_log_metadata[,grepl("ENSG", colnames(pGBM_log_metadata))] # Only genes!
colnames(pGBM_log) <- gsub("\\..*", "", colnames(pGBM_log))

# Load GO gene sets.
genesets <- read.csv(
    paste0(getwd(), "/data/c5.go.v7.5.1.symbols.gmt"),
    sep = "\t",
    header = FALSE,
    row.names = 1)

# Remove column that contains accession links & transpose datafrae.
genesets$V2 <- NULL 
genesets <- as.data.frame(t(genesets))

# Signature survival ------------------------------------------------------

signature_survival <- function(geneset.name) {

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
    colnames(signature.scores) <- c("scores","months","status")
    
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
        legend.labs = c(paste0("Above cutoff (n=", nrow(signature.scores[signature.scores$classification == "High",]), ")"), paste0("Below cutoff (n=", nrow(signature.scores[signature.scores$classification == "Low",]), ")")),
        xlab = "Time in months",
        title = geneset.name,
        risk.table = TRUE,
        legend.title = paste0("Cutoff = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))
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

    pdf(paste0(getwd(), "/results/TCGA-GBM/", geneset.name, "_survival_analysis.pdf"), width = 10, height = 6)
    print(surv.plot)
    dev.off()
}


## Calculate signature -----------------------------------------------------

signature_survival("GOCC_VESICLE_LUMEN") 
signature_survival("GOBP_EXTRACELLULAR_VESICLE_BIOGENESIS")
signature_survival("GOBP_EXOCYTOSIS")
signature_survival("GOBP_INTERCELLULAR_TRANSPORT")
signature_survival("GOBP_EXPORT_FROM_CELL")
signature_survival("GOCC_EXOCYTIC_VESICLE")
signature_survival("GOBP_REGULATION_OF_VESICLE_FUSION")
signature_survival("GOBP_EXOCYTIC_PROCESS")
signature_survival("GOBP_VESICLE_DOCKING")
signature_survival("GOCC_TRANSPORT_VESICLE")
signature_survival("GOCC_SECRETORY_VESICLE") 
signature_survival("GOCC_VESICLE_MEMBRANE") 
signature_survival("GOBP_ENDOCYTOSIS") 
signature_survival("GOBP_SECRETION")
signature_survival("GOBP_MEMBRANE_INVAGINATION") 
signature_survival("GOBP_VESICLE_ORGANIZATION")

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
    colnames(gene.exp) <- c("expression","months","status")
    
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
        legend.labs = c(paste0("Above cutoff (n=", nrow(gene.exp[gene.exp$classification == "High",]), ")"), paste0("Below cutoff (n=", nrow(gene.exp[gene.exp$classification == "Low",]), ")")),
        xlab = "Time in months",
        title = gene,
        risk.table = TRUE,
        legend.title = paste0("Cutoff = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))
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

    pdf(paste0(getwd(), "/results/TCGA-GBM/", gene, "_survival_analysis.pdf"), width = 10, height = 6)
    print(surv.plot)
    dev.off()
}

gene_survival("PRNP")
gene_survival("DSTN")
gene_survival("ANXA1")
gene_survival("RAB31")
gene_survival("SYPL1")