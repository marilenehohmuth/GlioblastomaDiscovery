# TCGA ANALYSES FOR PAPER ######################################################

# R version: 4.0.3
                                  # VERSIONS
library(edgeR)                    # 3.32.1
library(ggplot2)                  # 3.3.5 
library(dplyr)                    # 1.0.7 
library(ggpubr)                   # 0.4.0
library(biomaRt)                  # 2.46.3 
library(GSVA)                     # 1.38.2
library(Hmisc)                    #
library(RColorBrewer)             #

setwd("manuscript2022/")

# ANALYSIS OF OTHER SOLID TUMORS FROM TCGA ######################################

genesets <- read.csv("data/c5.go.v7.5.1.symbols.gmt",
                     header = FALSE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)
genesets.t <- t(genesets)
genesets.t <- as.data.frame(genesets.t)

traffic.genesets <- genesets.t %>% dplyr::select(GOCC_COATED_VESICLE,
                                                 GOCC_ENDOCYTIC_VESICLE,
                                                 GOCC_PHAGOCYTIC_VESICLE,
                                                 GOCC_SECRETORY_VESICLE,
                                                 GOBP_VESICLE_CYTOSKELETAL_TRAFFICKING,
                                                 GOBP_EXOCYTOSIS,
                                                 GOBP_ENDOCYTOSIS,
                                                 GOBP_SECRETION,
                                                 GOBP_INTRACELLULAR_TRANSPORT,
                                                 GOBP_INTRACELLULAR_PROTEIN_TRANSPORT,
                                                 GOBP_TRANSPORT_ALONG_MICROTUBULE,
                                                 GOCC_VESICLE_LUMEN,
                                                 GOCC_VESICLE_MEMBRANE,
                                                 GOBP_VESICLE_ORGANIZATION,
                                                 GOBP_VESICLE_TARGETING,
                                                 GOBP_VESICLE_DOCKING,
                                                 GOCC_CLATHRIN_COATED_VESICLE,
                                                 GOCC_COPI_COATED_VESICLE,
                                                 GOBP_VESICLE_TETHERING,
                                                 GOBP_INTERCELLULAR_TRANSPORT,
                                                 GOBP_REGULATION_OF_VESICLE_FUSION)

genesets.list <- list()
for (i in colnames(traffic.genesets)) {
  genesets.list[[tolower(i)]] <- traffic.genesets[,i][traffic.genesets[,i] != ""]
}

names(genesets.list) <- gsub(x = names(genesets.list),
                             pattern = c("gobp_"),
                             replacement = "")
names(genesets.list) <- gsub(x = names(genesets.list),
                             pattern = c("gocc_"),
                             replacement = "")
names(genesets.list) <- gsub(x = names(genesets.list),
                             pattern = c("_"),
                             replacement = " ")

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)

########## BREAST CANCER ##########

########## Step 1: Data retrieval ##########

BRCA_count_data <- readRDS("data/GDCdata/BRCA/TCGA-BRCA_1222_samples_raw_count_data.RDS")
                            
BRCA_metadata <- readRDS("data/GDCdata/BRCA/TCGA-BRCA_1222_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(BRCA_count_data),
                 mart = ensembl)
idx <- match(rownames(BRCA_count_data), genemap$ensembl_gene_id)
BRCA_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(BRCA_count_data) <- make.names(BRCA_count_data$hgnc_symbol, unique = TRUE)
BRCA_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
BRCA_primary_tumor_metadata <- BRCA_metadata[BRCA_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
BRCA_primary_tumor_count_data <- BRCA_count_data[,colnames(BRCA_count_data) %in% 
                                                 rownames(BRCA_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pBRCA_cpm <- cpm(BRCA_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pBRCA_cpm1 <- pBRCA_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pBRCA_log <- as.data.frame(log10(pBRCA_cpm1))

########## Step 6: Gene set variation analysis & correlation assessment ##########

BRCA.df <- data.frame()
for (stage in unique(BRCA_primary_tumor_metadata$ajcc_pathologic_stage[is.na(BRCA_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  print(stage)
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(BRCA_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame BRCAording to the given stage.
  stage.metadata <- BRCA_primary_tumor_metadata[which(BRCA_primary_tumor_metadata$ajcc_pathologic_stage == stage),]

  stage.norm.counts <- pBRCA_log[,colnames(pBRCA_log) %in% rownames(stage.metadata)]

  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)

  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))

  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("BRCA",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  BRCA.df <- rbind(BRCA.df, corr.prnp)
}                            

pBRCA_cpm <- NULL
pBRCA_cpm1 <- NULL
pBRCA_log <- NULL
BRCA_count_data <- NULL
BRCA_metadata <- NULL
BRCA_primary_tumor_count_data <- NULL
BRCA_primary_tumor_metadata <- NULL

##### COLON CANCER ######

########## Step 1: Data retrieval ##########

COAD_count_data <- readRDS("data/GDCdata/COAD/TCGA-COAD_521_samples_raw_count_data.RDS")
                            
COAD_metadata <- readRDS("data/GDCdata/COAD/TCGA-COAD_521_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(COAD_count_data),
                 mart = ensembl)
idx <- match(rownames(COAD_count_data), genemap$ensembl_gene_id)
COAD_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(COAD_count_data) <- make.names(COAD_count_data$hgnc_symbol, unique = TRUE)
COAD_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
COAD_primary_tumor_metadata <- COAD_metadata[COAD_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
COAD_primary_tumor_count_data <- COAD_count_data[,colnames(COAD_count_data) %in% 
                                                   rownames(COAD_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pCOAD_cpm <- cpm(COAD_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pCOAD_cpm1 <- pCOAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pCOAD_log <- as.data.frame(log10(pCOAD_cpm1))

########## Step 5: Gene set variation analysis & correlation assessment ##########

COAD.df <- data.frame()
for (stage in unique(COAD_primary_tumor_metadata$ajcc_pathologic_stage[is.na(COAD_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(COAD_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame COADording to the given stage.
  stage.metadata <- COAD_primary_tumor_metadata[which(COAD_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pCOAD_log[,colnames(pCOAD_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("COAD",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  COAD.df <- rbind(COAD.df, corr.prnp)
}        
pCOAD_cpm <- NULL
pCOAD_cpm1 <- NULL
pCOAD_log <- NULL
COAD_count_data <- NULL
COAD_metadata <- NULL
COAD_primary_tumor_count_data <- NULL
COAD_primary_tumor_metadata <- NULL

##### RECTUM CANCER ######

########## Step 1: Data retrieval ##########

READ_count_data <- readRDS("data/GDCdata/READ/TCGA-READ_177_samples_raw_count_data.RDS")
                            
READ_metadata <- readRDS("data/GDCdata/READ/TCGA-READ_177_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(READ_count_data),
                 mart = ensembl)
idx <- match(rownames(READ_count_data), genemap$ensembl_gene_id)
READ_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(READ_count_data) <- make.names(READ_count_data$hgnc_symbol, unique = TRUE)
READ_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
READ_primary_tumor_metadata <- READ_metadata[READ_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
READ_primary_tumor_count_data <- READ_count_data[,colnames(READ_count_data) %in% 
                                                   rownames(READ_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pREAD_cpm <- cpm(READ_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pREAD_cpm1 <- pREAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pREAD_log <- as.data.frame(log10(pREAD_cpm1))

print("Normalized primary tumor count data.")

########## Step 5: Gene set variation analysis & correlation assessment ##########

READ.df <- data.frame()
for (stage in unique(READ_primary_tumor_metadata$ajcc_pathologic_stage[is.na(READ_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(READ_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame READording to the given stage.
  stage.metadata <- READ_primary_tumor_metadata[which(READ_primary_tumor_metadata$ajcc_pathologic_stage == stage),]

  stage.norm.counts <- pREAD_log[,colnames(pREAD_log) %in% rownames(stage.metadata)]

  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)

  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))

  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("READ",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  READ.df <- rbind(READ.df, corr.prnp)
}                
pREAD_cpm <- NULL
pREAD_cpm1 <- NULL
pREAD_log <- NULL
READ_count_data <- NULL
READ_metadata <- NULL
READ_primary_tumor_count_data <- NULL
READ_primary_tumor_metadata <- NULL

##### PROSTATE CANCER ######

########## Step 1: Data retrieval ##########

PRAD_count_data <- readRDS("data/GDCdata/PRAD/TCGA-PRAD_551_samples_raw_count_data.RDS")
                            
PRAD_metadata <- readRDS("data/GDCdata/PRAD/TCGA-PRAD_551_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(PRAD_count_data),
                 mart = ensembl)
idx <- match(rownames(PRAD_count_data), genemap$ensembl_gene_id)
PRAD_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(PRAD_count_data) <- make.names(PRAD_count_data$hgnc_symbol, unique = TRUE)
PRAD_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
PRAD_primary_tumor_metadata <- PRAD_metadata[PRAD_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
PRAD_primary_tumor_count_data <- PRAD_count_data[,colnames(PRAD_count_data) %in% 
                                                   rownames(PRAD_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pPRAD_cpm <- cpm(PRAD_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pPRAD_cpm1 <- pPRAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pPRAD_log <- as.data.frame(log10(pPRAD_cpm1))

print("Normalized primary tumor count data.")

########## Step 5: Gene set variation analysis ##########

PRAD.gsva <- gsva(as.matrix(pPRAD_log),
                  genesets.list)

PRAD.prnp <- as.data.frame(t(pPRAD_log["PRNP",]))

PRAD.gsva.prnp <- cbind(t(PRAD.gsva), PRAD.prnp)
PRAD.gsva.prnp.metadata <- cbind(PRAD.gsva.prnp, PRAD_primary_tumor_metadata)

corr <- rcorr(as.matrix(PRAD.gsva.prnp.metadata[,1:22]),
              type = "spearman")
PRAD.df <- cbind(corr$r["PRNP",], 
                 corr$P["PRNP",])
colnames(PRAD.df) <- c("R", "P")
PRAD.df <- PRAD.df[-c(22),]
PRAD.df <- cbind(PRAD.df, 
                 Cancer.type = as.vector(rep("PRAD",  each = 21)), 
                 Process = rownames(corr.prnp),
                 Stage = as.vector(rep("Not reported",  each = 21)))
 
PRAD.df <- as.data.frame(PRAD.df)


pPRAD_cpm <- NULL
pPRAD_cpm1 <- NULL
pPRAD_log <- NULL
PRAD_count_data <- NULL
PRAD_metadata <- NULL
PRAD_primary_tumor_count_data <- NULL
PRAD_primary_tumor_metadata <- NULL

##### LUNG CANCER ######

########## Step 1: Data retrieval ##########

LUAD_count_data <- readRDS("data/GDCdata/LUAD/TCGA-LUAD_594_samples_raw_count_data.RDS")
                            
LUAD_metadata <- readRDS("data/GDCdata/LUAD/TCGA-LUAD_594_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(LUAD_count_data),
                 mart = ensembl)
idx <- match(rownames(LUAD_count_data), genemap$ensembl_gene_id)
LUAD_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(LUAD_count_data) <- make.names(LUAD_count_data$hgnc_symbol, unique = TRUE)
LUAD_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
LUAD_primary_tumor_metadata <- LUAD_metadata[LUAD_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
LUAD_primary_tumor_count_data <- LUAD_count_data[,colnames(LUAD_count_data) %in% 
                                                   rownames(LUAD_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pLUAD_cpm <- cpm(LUAD_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pLUAD_cpm1 <- pLUAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pLUAD_log <- as.data.frame(log10(pLUAD_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########

LUAD.df <- data.frame()
for (stage in unique(LUAD_primary_tumor_metadata$ajcc_pathologic_stage[is.na(LUAD_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(LUAD_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame LUADording to the given stage.
  stage.metadata <- LUAD_primary_tumor_metadata[which(LUAD_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pLUAD_log[,colnames(pLUAD_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("LUAD",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  LUAD.df <- rbind(LUAD.df, corr.prnp)
}        

pLUAD_cpm <- NULL
pLUAD_cpm1 <- NULL
pLUAD_log <- NULL
LUAD_count_data <- NULL
LUAD_metadata <- NULL
LUAD_primary_tumor_count_data <- NULL
LUAD_primary_tumor_metadata <- NULL

##### RENAL CANCER ######

########## Step 1: Data retrieval ##########

KIRC_count_data <- readRDS("data/GDCdata/KIRC/TCGA-KIRC_611_samples_raw_count_data.RDS")
                            
                            
KIRC_metadata <- readRDS("data/GDCdata/KIRC/TCGA-KIRC_611_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(KIRC_count_data),
                 mart = ensembl)
idx <- match(rownames(KIRC_count_data), genemap$ensembl_gene_id)
KIRC_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(KIRC_count_data) <- make.names(KIRC_count_data$hgnc_symbol, unique = TRUE)
KIRC_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
KIRC_primary_tumor_metadata <- KIRC_metadata[KIRC_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
KIRC_primary_tumor_count_data <- KIRC_count_data[,colnames(KIRC_count_data) %in% 
                                                   rownames(KIRC_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pKIRC_cpm <- cpm(KIRC_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pKIRC_cpm1 <- pKIRC_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pKIRC_log <- as.data.frame(log10(pKIRC_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########

KIRC.df <- data.frame()
for (stage in unique(KIRC_primary_tumor_metadata$ajcc_pathologic_stage[is.na(KIRC_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(KIRC_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame KIRCording to the given stage.
  stage.metadata <- KIRC_primary_tumor_metadata[which(KIRC_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pKIRC_log[,colnames(pKIRC_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("KIRC",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  KIRC.df <- rbind(KIRC.df, corr.prnp)
}            

pKIRC_cpm <- NULL
pKIRC_cpm1 <- NULL
pKIRC_log <- NULL
KIRC_count_data <- NULL
KIRC_metadata <- NULL
KIRC_primary_tumor_count_data <- NULL
KIRC_primary_tumor_metadata <- NULL

##### GASTRIC CANCER ######

########## Step 1: Data retrieval ##########

STAD_count_data <- readRDS("data/GDCdata/STAD/TCGA-STAD_407_samples_raw_count_data.RDS")
                            
STAD_metadata <- readRDS("data/GDCdata/STAD/TCGA-STAD_407_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",
                     mart = ensembl)
genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(STAD_count_data),
                 mart = ensembl)
idx <- match(rownames(STAD_count_data), genemap$ensembl_gene_id)
STAD_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(STAD_count_data) <- make.names(STAD_count_data$hgnc_symbol, unique = TRUE)
STAD_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
STAD_primary_tumor_metadata <- STAD_metadata[STAD_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
STAD_primary_tumor_count_data <- STAD_count_data[,colnames(STAD_count_data) %in% 
                                                   rownames(STAD_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pSTAD_cpm <- cpm(STAD_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pSTAD_cpm1 <- pSTAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pSTAD_log <- as.data.frame(log10(pSTAD_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########

STAD.df <- data.frame()
for (stage in unique(STAD_primary_tumor_metadata$ajcc_pathologic_stage[is.na(STAD_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(STAD_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame STADording to the given stage.
  stage.metadata <- STAD_primary_tumor_metadata[which(STAD_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pSTAD_log[,colnames(pSTAD_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("STAD",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  STAD.df <- rbind(STAD.df, corr.prnp)
}                            


pSTAD_cpm <- NULL
pSTAD_cpm1 <- NULL
pSTAD_log <- NULL
STAD_count_data <- NULL
STAD_metadata <- NULL
STAD_primary_tumor_count_data <- NULL
STAD_primary_tumor_metadata <- NULL

##### HEAD AND NECK CANCER ######

########## Step 1: Data retrieval ##########

HNSC_count_data <- readRDS("data/GDCdata/HNSC/TCGA-HNSC_546_samples_raw_count_data.RDS")
                            
                            
HNSC_metadata <- readRDS("data/GDCdata/HNSC/TCGA-HNSC_546_samples_metadata.RDS")
                          

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(HNSC_count_data),
                 mart = ensembl)
idx <- match(rownames(HNSC_count_data), genemap$ensembl_gene_id)
HNSC_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(HNSC_count_data) <- make.names(HNSC_count_data$hgnc_symbol, unique = TRUE)
HNSC_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
HNSC_primary_tumor_metadata <- HNSC_metadata[HNSC_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
HNSC_primary_tumor_count_data <- HNSC_count_data[,colnames(HNSC_count_data) %in% 
                                                   rownames(HNSC_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pHNSC_cpm <- cpm(HNSC_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pHNSC_cpm1 <- pHNSC_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pHNSC_log <- as.data.frame(log10(pHNSC_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


HNSC.df <- data.frame()
for (stage in unique(HNSC_primary_tumor_metadata$ajcc_pathologic_stage[is.na(HNSC_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(HNSC_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame HNSCording to the given stage.
  stage.metadata <- HNSC_primary_tumor_metadata[which(HNSC_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pHNSC_log[,colnames(pHNSC_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("HNSC",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  HNSC.df <- rbind(HNSC.df, corr.prnp)
}                            

pHNSC_cpm <- NULL
pHNSC_cpm1 <- NULL
pHNSC_log <- NULL
HNSC_count_data <- NULL
HNSC_metadata <- NULL
HNSC_primary_tumor_count_data <- NULL
HNSC_primary_tumor_metadata <- NULL

##### OVARIAN CANCER ######

########## Step 1: Data retrieval ##########

OV_count_data <- readRDS("data/GDCdata/OV/TCGA-OV_379_samples_raw_count_data.RDS")
                            
                            
OV_metadata <- readRDS("data/GDCdata/OV/TCGA-OV_379_samples_metadata.RDS")

########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(OV_count_data),
                 mart = ensembl)
idx <- match(rownames(OV_count_data), genemap$ensembl_gene_id)
OV_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(OV_count_data) <- make.names(OV_count_data$hgnc_symbol, unique = TRUE)
OV_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
OV_primary_tumor_metadata <- OV_metadata[OV_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
OV_primary_tumor_count_data <- OV_count_data[,colnames(OV_count_data) %in% 
                                                   rownames(OV_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pOV_cpm <- cpm(OV_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pOV_cpm1 <- pOV_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pOV_log <- as.data.frame(log10(pOV_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


OV.df <- data.frame()
for (stage in unique(OV_primary_tumor_metadata$figo_stage[is.na(OV_primary_tumor_metadata$figo_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(OV_primary_tumor_metadata$figo_stage)[stage] <= 4) next
  
  # Subset data.frame OVording to the given stage.
  stage.metadata <- OV_primary_tumor_metadata[which(OV_primary_tumor_metadata$figo_stage == stage),]
  print("ok")
  stage.norm.counts <- pOV_log[,colnames(pOV_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("OV",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  OV.df <- rbind(OV.df, corr.prnp)
}                            

pOV_cpm <- NULL
pOV_cpm1 <- NULL
pOV_log <- NULL
OV_count_data <- NULL
OV_metadata <- NULL
OV_primary_tumor_count_data <- NULL
OV_primary_tumor_metadata <- NULL

##### PANCREATIC CANCER ######

########## Step 1: Data retrieval ##########

PAAD_count_data <- readRDS("data/GDCdata/PAAD/TCGA-PAAD_182_samples_raw_count_data.RDS")
                            
                            
PAAD_metadata <- readRDS("data/GDCdata/PAAD/TCGA-PAAD_182_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(PAAD_count_data),
                 mart = ensembl)
idx <- match(rownames(PAAD_count_data), genemap$ensembl_gene_id)
PAAD_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(PAAD_count_data) <- make.names(PAAD_count_data$hgnc_symbol, unique = TRUE)
PAAD_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
PAAD_primary_tumor_metadata <- PAAD_metadata[PAAD_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
PAAD_primary_tumor_count_data <- PAAD_count_data[,colnames(PAAD_count_data) %in% 
                                                   rownames(PAAD_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pPAAD_cpm <- cpm(PAAD_primary_tumor_count_data, 
               log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pPAAD_cpm1 <- pPAAD_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pPAAD_log <- as.data.frame(log10(pPAAD_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########

PAAD.df <- data.frame()
for (stage in unique(PAAD_primary_tumor_metadata$ajcc_pathologic_stage[is.na(PAAD_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(PAAD_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame PAADording to the given stage.
  stage.metadata <- PAAD_primary_tumor_metadata[which(PAAD_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pPAAD_log[,colnames(pPAAD_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("PAAD",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  PAAD.df <- rbind(PAAD.df, corr.prnp)
}           

pPAAD_cpm <- NULL
pPAAD_cpm1 <- NULL
pPAAD_log <- NULL
PAAD_count_data <- NULL
PAAD_metadata <- NULL
PAAD_primary_tumor_count_data <- NULL
PAAD_primary_tumor_metadata <- NULL

##### UTERINE CANCER ######

########## Step 1: Data retrieval ##########

UCEC_count_data <- readRDS("data/GDCdata/UCEC/TCGA-UCEC_587_samples_raw_count_data.RDS")
                            
                            
UCEC_metadata <- readRDS("data/GDCdata/UCEC/TCGA-UCEC_587_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(UCEC_count_data),
                 mart = ensembl)
idx <- match(rownames(UCEC_count_data), genemap$ensembl_gene_id)
UCEC_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(UCEC_count_data) <- make.names(UCEC_count_data$hgnc_symbol, unique = TRUE)
UCEC_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
UCEC_primary_tumor_metadata <- UCEC_metadata[UCEC_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
UCEC_primary_tumor_count_data <- UCEC_count_data[,colnames(UCEC_count_data) %in% 
                                                   rownames(UCEC_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pUCEC_cpm <- cpm(UCEC_primary_tumor_count_data, 
               log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pUCEC_cpm1 <- pUCEC_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pUCEC_log <- as.data.frame(log10(pUCEC_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


UCEC.df <- data.frame()
for (stage in unique(UCEC_primary_tumor_metadata$figo_stage[is.na(UCEC_primary_tumor_metadata$figo_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(UCEC_primary_tumor_metadata$figo_stage)[stage] <= 4) next
  
  # Subset data.frame UCECording to the given stage.
  stage.metadata <- UCEC_primary_tumor_metadata[which(UCEC_primary_tumor_metadata$figo_stage == stage),]
  print("ok")
  stage.norm.counts <- pUCEC_log[,colnames(pUCEC_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("UCEC",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  UCEC.df <- rbind(UCEC.df, corr.prnp)
}                            

pUCEC_log <- NULL
UCEC_count_data <- NULL
UCEC_metadata <- NULL

##### BLADDER CANCER ######

########## Step 1: Data retrieval ##########

BLCA_count_data <- readRDS("data/GDCdata/BLCA/TCGA-BLCA_433_samples_raw_count_data.RDS")
                            
                            
BLCA_metadata <- readRDS("data/GDCdata/BLCA/TCGA-BLCA_433_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(BLCA_count_data),
                 mart = ensembl)
idx <- match(rownames(BLCA_count_data), genemap$ensembl_gene_id)
BLCA_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(BLCA_count_data) <- make.names(BLCA_count_data$hgnc_symbol, unique = TRUE)
BLCA_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
BLCA_primary_tumor_metadata <- BLCA_metadata[BLCA_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
BLCA_primary_tumor_count_data <- BLCA_count_data[,colnames(BLCA_count_data) %in% 
                                                   rownames(BLCA_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pBLCA_cpm <- cpm(BLCA_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pBLCA_cpm1 <- pBLCA_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pBLCA_log <- as.data.frame(log10(pBLCA_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


BLCA.df <- data.frame()
for (stage in unique(BLCA_primary_tumor_metadata$ajcc_pathologic_stage[is.na(BLCA_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(BLCA_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame BLCAording to the given stage.
  stage.metadata <- BLCA_primary_tumor_metadata[which(BLCA_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pBLCA_log[,colnames(pBLCA_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("BLCA",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  BLCA.df <- rbind(BLCA.df, corr.prnp)
}        

pBLCA_log <- NULL
BLCA_count_data <- NULL
BLCA_metadata <- NULL

##### GALLBLADDER CANCER ######

########## Step 1: Data retrieval ##########

CHOL_count_data <- readRDS("data/GDCdata/CHOL/TCGA-CHOL_45_samples_raw_count_data.RDS")
                            
                            
CHOL_metadata <- readRDS("data/GDCdata/CHOL/TCGA-CHOL_45_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(CHOL_count_data),
                 mart = ensembl)
idx <- match(rownames(CHOL_count_data), genemap$ensembl_gene_id)
CHOL_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(CHOL_count_data) <- make.names(CHOL_count_data$hgnc_symbol, unique = TRUE)
CHOL_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
CHOL_primary_tumor_metadata <- CHOL_metadata[CHOL_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
CHOL_primary_tumor_count_data <- CHOL_count_data[,colnames(CHOL_count_data) %in% 
                                                   rownames(CHOL_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pCHOL_cpm <- cpm(CHOL_primary_tumor_count_data, 
                 log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pCHOL_cpm1 <- pCHOL_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pCHOL_log <- as.data.frame(log10(pCHOL_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


CHOL.df <- data.frame()
for (stage in unique(CHOL_primary_tumor_metadata$ajcc_pathologic_stage[is.na(CHOL_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(CHOL_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame CHOLording to the given stage.
  stage.metadata <- CHOL_primary_tumor_metadata[which(CHOL_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pCHOL_log[,colnames(pCHOL_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("CHOL",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  CHOL.df <- rbind(CHOL.df, corr.prnp)
} 

pCHOL_log <- NULL
CHOL_count_data <- NULL
CHOL_metadata <- NULL

##### ADRENAL GLAND CANCER ######

########## Step 1: Data retrieval ##########

ACC_count_data <- readRDS("data/GDCdata/ACC/TCGA-ACC_79_samples_raw_count_data.RDS")
                            
                            
ACC_metadata <- readRDS("data/GDCdata/ACC/TCGA-ACC_79_samples_metadata.RDS")
                          
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(ACC_count_data),
                 mart = ensembl)
idx <- match(rownames(ACC_count_data), genemap$ensembl_gene_id)
ACC_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(ACC_count_data) <- make.names(ACC_count_data$hgnc_symbol, unique = TRUE)
ACC_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
ACC_primary_tumor_metadata <- ACC_metadata[ACC_metadata$definition == "Primary solid Tumor",]
  
# Filtering raw count data.
ACC_primary_tumor_count_data <- ACC_count_data[,colnames(ACC_count_data) %in% 
                                                 rownames(ACC_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pACC_cpm <- cpm(ACC_primary_tumor_count_data, 
                log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pACC_cpm1 <- pACC_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pACC_log <- as.data.frame(log10(pACC_cpm1))

########## Step 6: GSVA ##########

ACC.df <- data.frame()
for (stage in unique(ACC_primary_tumor_metadata$ajcc_pathologic_stage[is.na(ACC_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(ACC_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame according to the given stage.
  stage.metadata <- ACC_primary_tumor_metadata[which(ACC_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pACC_log[,colnames(pACC_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("ACC",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  ACC.df <- rbind(ACC.df, corr.prnp)
}                            
            
pACC_log <- NULL
ACC_count_data <- NULL
ACC_metadata <- NULL

##### MELANOMA ######

########## Step 1: Data retrieval ##########

SKCM_count_data <- readRDS("data/GDCdata/SKCM/TCGA-SKCM_472_samples_raw_count_data.RDS")
                            
SKCM_metadata <- readRDS("data/GDCdata/SKCM/TCGA-SKCM_472_samples_metadata.RDS")
    
########## Step 2: Converting ENSEMBL IDs to gene symbols ##########

genemap <- getBM(filters = c("ensembl_gene_id"),
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = rownames(SKCM_count_data),
                 mart = ensembl)
idx <- match(rownames(SKCM_count_data), genemap$ensembl_gene_id)
SKCM_count_data$hgnc_symbol <- genemap$hgnc_symbol[idx]

rownames(SKCM_count_data) <- make.names(SKCM_count_data$hgnc_symbol, unique = TRUE)
SKCM_count_data$hgnc_symbol <- NULL

########## Step 3: Selecting primary tumors ##########

# Filtering metadata.
SKCM_primary_tumor_metadata <- SKCM_metadata[SKCM_metadata$definition == "Primary solid Tumor",]

# Filtering raw count data.
SKCM_primary_tumor_count_data <- SKCM_count_data[,colnames(SKCM_count_data) %in% 
                                                   rownames(SKCM_primary_tumor_metadata)]

########## Step 4: Normalizing primary tumor count data ##########

# Normalizing raw count data of primary tumors (CPM normalization).
pSKCM_cpm <- cpm(SKCM_primary_tumor_count_data, 
                log = FALSE)

# Summing 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
pSKCM_cpm1 <- pSKCM_cpm + 1

# Log-transforming CPM-normalized count data of primary tumors.
pSKCM_log <- as.data.frame(log10(pSKCM_cpm1))

print("Normalized primary tumor count data.")

########## Step 6: GSVA ##########


SKCM.df <- data.frame()
for (stage in unique(SKCM_primary_tumor_metadata$ajcc_pathologic_stage[is.na(SKCM_primary_tumor_metadata$ajcc_pathologic_stage) == FALSE])) {
  
  # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis)
  if(table(SKCM_primary_tumor_metadata$ajcc_pathologic_stage)[stage] <= 4) next
  
  # Subset data.frame SKCMording to the given stage.
  stage.metadata <- SKCM_primary_tumor_metadata[which(SKCM_primary_tumor_metadata$ajcc_pathologic_stage == stage),]
  print("ok")
  stage.norm.counts <- pSKCM_log[,colnames(pSKCM_log) %in% rownames(stage.metadata)]
  print("ok2")
  stage.gsva <- gsva(as.matrix(stage.norm.counts),
                     genesets.list)
  print("ok3")
  stage.gsva.prnp <- cbind(t(stage.gsva), as.data.frame(t(stage.norm.counts["PRNP",])))
  print("ok4")
  stage.gsva.prnp.metadata <- cbind(stage.gsva.prnp, stage.metadata)
  print("ok5")
  
  # Calculate Spearman correlations in the matrix
  corr <- rcorr(as.matrix(stage.gsva.prnp.metadata[,1:22]),
                type = "spearman")
  
  # Put correlation coefficients and p-values related to PRNP in a dataframe.  
  corr.prnp <- cbind(corr$r["PRNP",], 
                     corr$P["PRNP",])
  colnames(corr.prnp) <- c("R", "P")
  corr.prnp <- corr.prnp[-c(22),] # Exclude because this line is PRNP vs PRNP.
  
  s <- strsplit(stage, split = " ")
  corr.prnp <- cbind(corr.prnp, 
                     Cancer.type = as.vector(rep("SKCM",  each = 21)), 
                     Process = rownames(corr.prnp),
                     Stage = as.vector(rep(s[[1]][2],  each = 21)))
  SKCM.df <- rbind(SKCM.df, corr.prnp)
}                            

pSKCM_log <- NULL
SKCM_count_data <- NULL
SKCM_metadata <- NULL

##################

df <- rbind(BRCA.df,
            COAD.df,
            READ.df,
            PRAD.df,
            LUAD.df,
            KIRC.df,
            STAD.df,
            HNSC.df,
            OV.df,
            PAAD.df,
            UCEC.df,
            BLCA.df,
            CHOL.df,
            ACC.df,
            SKCM.df)
write.csv(df,
          file = "results/TCGA/other_tumors/TCGA_correlation_PRNPvsTrafficGenesets.csv")


df <- read.csv("results/TCGA/other_tumors/TCGA_correlation_PRNPvsTrafficGenesets.csv",
               header = TRUE,
               row.names = 1)

labs <- c("Breast", "Colon", "Rectum", "Prostate", "Lung", "Kidney",
               "Stomach", "Head and neck", "Ovary", "Pancreas", "Uterus",
               "Bladder", "Gallbladder", "Adrenal gland", "Skin")
names(labs) <- c("BRCA", "COAD", "READ",
                      "PRAD", "LUAD", "KIRC",
                      "STAD", "HNSC", "OV",
                      "PAAD", "UCEC", "BLCA",
                      "CHOL", "ACC", "SKCM")

pdf("results/TCGA/other_tumors/TCGA_correlation_PRNPvsTrafficGenesets.pdf",
    width = 6.5,
    height = 15)
ggplot(df,
       aes(y = Stage, x = Process)) +
  geom_point(aes(fill = as.numeric(R), 
                 size = -log10(as.numeric(P))),
             shape = 21,
             color = "black") +
  scale_fill_gradientn(limits = c(min(as.numeric(df$R)), max(as.numeric(df$R))),
                       colors = colorRampPalette(rev(brewer.pal(9, "RdBu")))(1000),
                       name = "Spearman correlation") +
  scale_size_continuous(limits = c(-log10(0.05),30),
                        name = expression(-log[10](p-value))) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1, 'cm'),
        axis.text.x = element_text(angle = 35,
                                   hjust = 1,
                                   vjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + 
  guides(fill = guide_colorbar(title.position = "top"),
         size = guide_legend(title.position = "top")) +
  facet_grid(rows = vars(Cancer.type), 
             scales = "free",
             space = "free", 
             labeller = labeller(Cancer.type = labs)) +
  xlab(NULL) +
  ylab(NULL)
dev.off()

sessionInfo()
