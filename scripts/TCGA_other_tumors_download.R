# DOWNLOAD OF OTHER TCGA PROJECTS ###############################################################
# @ This script downloads RNA-seq data of the following TCGA projects:                          # 
# @ BRCA, COAD, READ, PRAD, LUAD, KIRC, STAD, HNSC, OV, PAAD, UCEC, BLCA, CHOL, ACC, PCPG, SKCM #                   
#################################################################################################


##########################
#### Loading packages ####
##########################
                                  
library(TCGAbiolinks)             # 2.23.5
library(SummarizedExperiment)     # 1.20.0


#################################################
#### Define a function to download TCGA data ####
#################################################

download_tcga <- function(project_id) {

    # Create query.
    query <- GDCquery(
        project = paste0("TCGA-", project_id),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "HTSeq - Counts",
        experimental.strategy = "RNA-Seq"
    )
    # Download data.
    GDCdownload(
        query,
        method = "api",
        files.per.chunk = 10
    )
    data <- GDCprepare(query)

    # Get count data.
    count_data <- as.data.frame(assay(data))

    # Get metadata.
    metadata <- as.data.frame(colData(data))

    # Get number of samples.
    n_samples <- nrow(metadata)

    # Define prefix for output file names.  
    file_prefix <- paste0("data/other_tumors/", project_id, "/TCGA-", project_id, "_", n_samples, "_samples")
    
    # Save count data to output file.
    saveRDS(
        count_data,
        file = paste0(file_prefix, "_count_data.RDS")
    )

    # Save metadata to output file.
    saveRDS(
        metadata,
        file = paste0(file_prefix, "_metadata.RDS")
    )
}


##############################################
#### Download data from each TCGA project ####
##############################################

# Define TCGA project IDs.
project_list <- c(
    "BRCA", "COAD", "READ", "PRAD", "LUAD", "KIRC", "STAD", "HNSC", 
    "OV", "PAAD", "UCEC", "BLCA", "CHOL", "ACC", "PCPG", "SKCM"
)

# Download TCGA projects that were specified.
for(project in project_list) {
    download_tcga(project_id = project)
}