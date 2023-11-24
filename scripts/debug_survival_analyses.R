
# Load packages -----------------------------------------------------------

library(survminer)                # 0.4-9
library(ggplot2)
library(dplyr)
library(GSVA)                     # 1.38.1
library(maxstat)                  # 0.7-25
library(survival)                 # 3.5-7
library(org.Hs.eg.db)             # 3.17.0



# Load data ---------------------------------------------------------------

pGBM_log_metadata <- readRDS("results/TCGA-GBM/PrimaryGBMs_logData+Metadata.RDS")
genesets <- read.csv(
  "data/c5.go.v7.5.1.symbols.gmt",
  sep = "\t",
  header = FALSE,
  row.names = 1)
genesets <- as.data.frame(t(genesets))

pGBM_log_metadata <- pGBM_log_metadata[pGBM_log_metadata$paper_IDH.status == "WT",]



# Survival analysis -------------------------------------------------------


##  Fix names ---------------------------------------------------------------

## removing version identifier from ENSEMBL IDs.
ensg_idx <- grep("^ENSG", colnames(pGBM_log_metadata))
start_ensg <- ensg_idx[1]
gene_ids <- sub("\\..*", "", colnames(pGBM_log_metadata)[ensg_idx])
colnames(pGBM_log_metadata)[ensg_idx] <- gene_ids

## getting correspondence between ENSEMBL IDs and gene symbols.
correspondence <- select(
  org.Hs.eg.db, 
  keys = gene_ids, 
  columns = "SYMBOL", 
  keytype = "ENSEMBL"
)

## replace column names
for (i in seq_len(nrow(correspondence))) {
  names(pGBM_log_metadata)[names(pGBM_log_metadata) == correspondence$ENSEMBL[i]] <- correspondence$SYMBOL[i]
}


## Run analysis ------------------------------------------------------------

geneset.name <- 'GOCC_VESICLE_LUMEN'

## function
genes.list <- list(geneset.name = genesets[,geneset.name][genesets[,geneset.name] != ""])
## bug
signature.scores <- gsva(expr = Matrix::t(pGBM_log_metadata[, start_ensg:length(names(pGBM_log_metadata))]), 
                         gset.idx.list = genes.list)
signature.scores <- as.data.frame(t(signature.scores))
signature.scores <- cbind(signature.scores, pGBM_log_metadata %>% select(paper_Survival..months., paper_Vital.status..1.dead.))
colnames(signature.scores) <- c("scores","months","status")
cutoff.test <- maxstat.test(
  Surv(months, status) ~ scores, 
  data = signature.scores,
  smethod = "LogRank",
  pmethod = "HL"
)
signature.scores <- signature.scores %>% mutate(classification = case_when(
  scores <= as.numeric(cutoff.test[["estimate"]]) ~ "Low",
  scores > as.numeric(cutoff.test[["estimate"]]) ~ "High"
))
fit <- survfit(
  Surv(months, status) ~ classification, 
  data = signature.scores,
  type = "kaplan-meier"
)

dist <- ggplot(
  signature.scores,
  aes(x = "", y = scores)
) +
  geom_boxplot() +
  geom_jitter(aes(color = classification)) +
  theme_classic() +
  ylab("GSVA scores") +
  xlab("") +
  theme(axis.ticks.x = element_blank(), legend.position = "top") +
  scale_color_discrete(name = "Signature")

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
  legend.labs = c("Above cutoff", "Below cutoff"),
  xlab = "Time in months",
  title = geneset.name,
  legend.title = paste0("Cutoff = ", signif(as.numeric(cutoff.test[["estimate"]]), digits = 4))
) 
plot.list <- list(
  "Dist" = dist,
  "Surv" = ggpubr::ggpar(surv.plot, font.legend = list(size = 15), font.x = list(size = 15),
                         font.y = list(size = 15), font.xtickslab = list(size = 15), font.ytickslab = list(size = 15))
)
return(plot.list)

