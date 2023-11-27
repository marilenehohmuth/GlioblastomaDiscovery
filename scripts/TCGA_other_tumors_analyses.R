# ANALYSIS OF OTHER TCGA PROJECTS ###############################################################
# @ This script does the analysis of RNA-seq data of the following TCGA projects:               # 
# @ BRCA, COAD, READ, PRAD, LUAD, KIRC, STAD, HNSC, OV, PAAD, UCEC, BLCA, CHOL, ACC, PCPG, SKCM #                   
#################################################################################################


##########################
#### Loading packages ####
##########################

library(edgeR)                    # 3.32.1
library(ggplot2)                  # 3.3.5 
library(dplyr)                    # 1.0.7 
library(ggpubr)                   # 0.4.0
library(GSVA)                     # 1.38.2
library(Hmisc)                    # 5.1_1
library(RColorBrewer)             # 1.1_3
library(org.Hs.eg.db)             # 3.17.0


########################################################################
#### Defining a function to process the data of other TCGA projects ####
########################################################################

process_tcga <- function(
  project_id, # TCGA project ID.
  genesets # List of gene sets to be considered in gene set variation analysis (GSVA).
) {
  
  # Read count data.
  count_data <- readRDS(paste0(getwd(), "/data/other_tumors/", project_id, "/TCGA-", project_id, "_count_data.RDS"))
  count_data$ensembl <- gsub("\\..*", "", rownames(count_data))

  # Read metadata.
  metadata <- readRDS(paste0(getwd(), "/data/other_tumors/", project_id, "/TCGA-", project_id, "_metadata.RDS"))

  # Get correspondence between ENSEMBL IDs and gene symbols.
  correspondence <- AnnotationDbi::select(
    org.Hs.eg.db, 
    keys = count_data$ensembl, 
    columns = "SYMBOL", 
    keytype = "ENSEMBL"
  )

  # Add gene symbols to count_data.
  count_data$hgnc_symbol <- correspondence$SYMBOL[match(count_data$ensembl, correspondence$ENSEMBL)]

  # Add gene symbols as row names, allowing duplicate row names by appending a suffix to repeated gene symbols.
  rownames(count_data) <- make.names(count_data$hgnc_symbol, unique = TRUE)
  
  # Remove column that has gene symbols, as these have been added as row names.
  count_data$hgnc_symbol <- NULL
  count_data$ensembl <- NULL

  # Keep only primary tumours in the metadata table.
  metadata <- metadata %>% dplyr::filter(definition == "Primary solid Tumor")

  # Keep only primary tumours in the count data table.
  count_data <- count_data %>% dplyr::select(rownames(metadata))
  
  # CPM-normalize count data.
  norm_count_data <- cpm(count_data, log = FALSE)

  # Sum 1 to all CPM-normalized counts in order to prevent issues with log-transformation.
  norm_count_data <- norm_count_data + 1

  # Log-transforming CPM-normalized count data.
  norm_count_data <- as.data.frame(log10(norm_count_data))

  # Adjust the name of the column that contains tumor pathological stages.
  if("ajcc_pathologic_stage" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "ajcc_pathologic_stage"] <- "tumor_stage"
  } else if("figo_stage" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "figo_stage"] <- "tumor_stage"
  } else if(nrow(metadata[grep("stage", colnames(metadata)),]) == 0) {
    metadata$tumor_stage <- NA
  }

  # Iterate over each pathological stage available for the given tumor...
  project_df <- data.frame()
  for(stage in metadata$tumor_stage[!is.na(metadata$tumor_stage)] %>% unique) {

    # Skip stage if there are only 4 samples or less (this number is not enough for correlation analysis).
    if(table(metadata$tumor_stage)[stage] <= 4) next

    # Subset metadata to the given stage.
    stage_metadata <- metadata %>% filter(tumor_stage == stage)
    
    # Subset normalized count data to the given stage.
    stage_norm_count_data <- norm_count_data %>% dplyr::select(rownames(stage_metadata))
    
    # Calculate GSVA scores for the given gene sets.
    stage_gsva <- gsva(as.matrix(stage_norm_count_data), genesets)

    # Add PRNP to GSVA table.
    stage_gsva <- cbind(t(stage_gsva), as.data.frame(t(stage_norm_count_data["PRNP",])))

    # Add metadata to GSVA table.
    stage_gsva <- cbind(stage_gsva, stage_metadata)

    # Number of columns to be included in the correlation analysis below (number of gene sets + PRNP column).
    n_columns <- length(genesets) + 1

    # Calculate Spearman correlations in the matrix
    corr <- rcorr(as.matrix(stage_gsva[,1:n_columns]), type = "spearman")

    # Put correlation coefficients and p-values related to PRNP in a dataframe.  
    corr_prnp <- cbind(corr$r["PRNP",], corr$P["PRNP",])
    colnames(corr_prnp) <- c("correlation", "pvalue")
    corr_prnp <- corr_prnp[rownames(corr_prnp) != "PRNP",] # Exclude because this line is PRNP vs PRNP.

    # Add more information to the dataframe.
    s <- strsplit(stage, split = " ")
    corr_prnp <- cbind(
      corr_prnp, 
      tumor = as.vector(rep(project_id, times = nrow(corr_prnp))), 
      stage = as.vector(rep(s[[1]][2], times = nrow(corr_prnp))),
      process = rownames(corr_prnp)
    )


    # Add stage results to the project-wide dataframe.
    project_df <- rbind(project_df, corr_prnp)
  }

  # Return project-wide dataframe containing the correlation between PRNP and the given gene sets.
  return(project_df)
}                            


###########################
#### Prepare gene sets ####
###########################

# Read gene sets file.
genesets <- read.csv(
  paste0(getwd(), "/data/c5.go.v7.5.1.symbols.gmt"),
  header = FALSE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)
# Transpose dataframe to have gene set names as column names.
genesets <- as.data.frame(t(genesets))

# Select gene sets related to traffic and vesicles.
traffic_genesets <- genesets %>% dplyr::select(
  GOCC_COATED_VESICLE, GOCC_ENDOCYTIC_VESICLE, GOCC_PHAGOCYTIC_VESICLE, GOCC_SECRETORY_VESICLE,
  GOBP_VESICLE_CYTOSKELETAL_TRAFFICKING, GOBP_EXOCYTOSIS, GOBP_ENDOCYTOSIS, GOBP_SECRETION,
  GOBP_INTRACELLULAR_TRANSPORT, GOBP_INTRACELLULAR_PROTEIN_TRANSPORT, GOBP_TRANSPORT_ALONG_MICROTUBULE,
  GOCC_VESICLE_LUMEN, GOCC_VESICLE_MEMBRANE, GOBP_VESICLE_ORGANIZATION, GOBP_VESICLE_TARGETING,
  GOBP_VESICLE_DOCKING, GOCC_CLATHRIN_COATED_VESICLE, GOCC_COPI_COATED_VESICLE, GOBP_VESICLE_TETHERING,
  GOBP_INTERCELLULAR_TRANSPORT, GOBP_REGULATION_OF_VESICLE_FUSION
)

# Change gene set names from upper to lower case.
genesets_list <- list()
for (column in colnames(traffic_genesets)) {
  genesets_list[[tolower(column)]] <- traffic_genesets[,column][traffic_genesets[,column] != ""]
}
# Further refine gene set names.
names(genesets_list) <- gsub(x = names(genesets_list), pattern = c("gobp_"), replacement = "")
names(genesets_list) <- gsub(x = names(genesets_list), pattern = c("gocc_"), replacement = "")
names(genesets_list) <- gsub(x = names(genesets_list), pattern = c("_"), replacement = " ")


####################################################
#### Processing the data of other TCGA projects ####
####################################################

# Define TCGA project IDs.
project_list <- c(
    "BRCA", "COAD", "READ", "PRAD", "LUAD", "KIRC", "STAD", "HNSC", "ESCA", "LIHC",
    "OV", "PAAD", "UCEC", "BLCA", "CHOL", "ACC", "PCPG", "SKCM", "THCA", "CESC"
)

# Process each TCGA project...
corr_all_projects <- data.frame()
for(project in project_list) {
    corr_project <- process_tcga(
      project_id = project,
      genesets = genesets_list
    )
    corr_all_projects <- rbind(corr_all_projects, corr_project)
}

# Save results to output file.
write.csv(corr_all_projects, file = paste0(getwd(), "/results/TCGA_other/TCGA_correlation_PRNPvsTrafficGenesets.csv"))


##########################
#### Plotting results ####
##########################

labs <- c(
  "Breast", "Colon", "Rectum", "Prostate", "Lung", "Kidney", "Stomach", "Head and neck", "Ovary", "Cervix", "Nerve/adrenal gland",
  "Pancreas", "Uterus", "Bladder", "Gallbladder", "Adrenal gland", "Skin", "Esophagus", "Liver", "Thyroid"
)
names(labs) <- c(
  "BRCA", "COAD", "READ", "PRAD", "LUAD", "KIRC", "STAD", "HNSC", "OV", "CESC", "PCPG",
  "PAAD", "UCEC", "BLCA", "CHOL", "ACC", "SKCM", "ESCA", "LIHC", "THCA")

corr_all_projects$stage[is.na(corr_all_projects$stage)] <- "Not reported"

#### @ FIGURE S8 (SUPPLEMENTAL) @ #### 
# Plotting the correlation between PRNP expression and traffic/vesicle signatures in
# various tumour types from TCGA.
pdf(paste0(getwd(), "/results/TCGA_other/TCGA_correlation_PRNPvsTrafficGenesets.pdf"),
    width = 6.5,
    height = 20)
ggplot(
  corr_all_projects[corr_all_projects$process != "PRNP",],
  aes(x = process, y = stage)
) +
  geom_point(
    aes(fill = as.numeric(correlation),  size = -log10(as.numeric(pvalue))), shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(min(as.numeric(corr_all_projects$correlation)), max(as.numeric(corr_all_projects$correlation))),
    colors = colorRampPalette(rev(brewer.pal(9, "RdBu")))(1000),
    name = "Spearman correlation"
  ) +
  scale_size_continuous(range = c(0.5,5), limits = c(-log10(0.05), 30), name = expression(-log[10](p-value))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(1, 'cm'),
    strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  ) + 
  guides(fill = guide_colorbar(title.position = "top"), size = guide_legend(title.position = "top")) +
  facet_grid(
    rows = vars(tumor), 
    scales = "free",
    space = "free", 
    labeller = labeller(tumor = labs)
  ) +
  xlab("Term") +
  ylab("Stage")
dev.off()