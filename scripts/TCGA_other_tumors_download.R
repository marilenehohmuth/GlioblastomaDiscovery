# TCGA ANALYSES FOR PAPER ######################################################

# R version: 4.0.3
                                  # VERSIONS
library(TCGAbiolinks)             # 2.23.5
library(SummarizedExperiment)     # 1.20.0

setwd("manuscript2022/")

########## BREAST CANCER ##########

query_BRCA <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      experimental.strategy = "RNA-Seq",
                      legacy = FALSE)
GDCdownload(query_BRCA,
            method = "api",
            files.per.chunk = 10)
data_BRCA <- GDCprepare(query_BRCA)

# Raw count data.
BRCA_count_data <- as.data.frame(assay(data_BRCA))
saveRDS(BRCA_count_data,
          file = "data/other_tumors/BRCA/TCGA-BRCA_1222_samples_raw_count_data.RDS")

# Metadata.
BRCA_metadata <- as.data.frame(colData(data_BRCA))
saveRDS(BRCA_metadata,
          file = "data/other_tumors/BRCA/TCGA-BRCA_1222_samples_metadata.RDS")

##### COLON CANCER ######

query_COAD <- GDCquery(project = "TCGA-COAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_COAD,
            method = "api",
            files.per.chunk = 10)
data_COAD <- GDCprepare(query_COAD)

# Raw count data.
COAD_count_data <- as.data.frame(assay(data_COAD))
saveRDS(COAD_count_data,
          file = "data/other_tumors/COAD/TCGA-COAD_521_samples_raw_count_data.RDS")

# Metadata.
COAD_metadata <- as.data.frame(colData(data_COAD))
saveRDS(COAD_metadata,
          file = "data/other_tumors/COAD/TCGA-COAD_521_samples_metadata.RDS")
  
##### RECTUM CANCER ######

query_READ <- GDCquery(project = "TCGA-READ",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_READ,
            method = "api",
            files.per.chunk = 10)
data_READ <- GDCprepare(query_READ)

# Raw count data.
READ_count_data <- as.data.frame(assay(data_READ))
saveRDS(READ_count_data,
          file = "data/other_tumors/READ/TCGA-READ_177_samples_raw_count_data.RDS")

# Metadata.
READ_metadata <- as.data.frame(colData(data_READ))
saveRDS(READ_metadata,
          file = "data/other_tumors/READ/TCGA-READ_177_samples_metadata.RDS")

##### PROSTATE CANCER ######

query_PRAD <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_PRAD,
            method = "api",
            files.per.chunk = 10)
data_PRAD <- GDCprepare(query_PRAD)

# Raw count data.
PRAD_count_data <- as.data.frame(assay(data_PRAD))
saveRDS(PRAD_count_data,
          file = "data/other_tumors/PRAD/TCGA-PRAD_551_samples_raw_count_data.RDS")

# Metadata.
PRAD_metadata <- as.data.frame(colData(data_PRAD))
saveRDS(PRAD_metadata,
          file = "data/other_tumors/PRAD/TCGA-PRAD_551_samples_metadata.RDS")

##### LUNG CANCER ######

query_LUAD <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_LUAD,
            method = "api",
            files.per.chunk = 10)
data_LUAD <- GDCprepare(query_LUAD)

# Raw count data.
LUAD_count_data <- as.data.frame(assay(data_LUAD))
saveRDS(LUAD_count_data,
          file = "data/other_tumors/LUAD/TCGA-LUAD_594_samples_raw_count_data.RDS")

# Metadata.
LUAD_metadata <- as.data.frame(colData(data_LUAD))
saveRDS(LUAD_metadata,
          file = "data/other_tumors/LUAD/TCGA-LUAD_594_samples_metadata.RDS")

##### RENAL CANCER ######

query_KIRC <- GDCquery(project = "TCGA-KIRC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_KIRC,
            method = "api",
            files.per.chunk = 10)
data_KIRC <- GDCprepare(query_KIRC)

# Raw count data.
KIRC_count_data <- as.data.frame(assay(data_KIRC))
saveRDS(KIRC_count_data,
          file = "data/other_tumors/KIRC/TCGA-KIRC_611_samples_raw_count_data.RDS")

# Metadata.
KIRC_metadata <- as.data.frame(colData(data_KIRC))
saveRDS(KIRC_metadata,
          file = "data/other_tumors/KIRC/TCGA-KIRC_611_samples_metadata.RDS")

##### GASTRIC CANCER ######

query_STAD <- GDCquery(project = "TCGA-STAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_STAD,
            method = "api",
            files.per.chunk = 10)
data_STAD <- GDCprepare(query_STAD)

# Raw count data.
STAD_count_data <- as.data.frame(assay(data_STAD))
saveRDS(STAD_count_data,
          file = "data/other_tumors/STAD/TCGA-STAD_407_samples_raw_count_data.RDS")

# Metadata.
STAD_metadata <- as.data.frame(colData(data_STAD))
saveRDS(STAD_metadata,
          file = "data/other_tumors/STAD/TCGA-STAD_407_samples_metadata.RDS")

##### HEAD AND NECK CANCER ######

query_HNSC <- GDCquery(project = "TCGA-HNSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_HNSC,
            method = "api",
            files.per.chunk = 10)
data_HNSC <- GDCprepare(query_HNSC)

# Raw count data.
HNSC_count_data <- as.data.frame(assay(data_HNSC))
saveRDS(HNSC_count_data,
          file = "data/other_tumors/HNSC/TCGA-HNSC_546_samples_raw_count_data.RDS")

# Metadata.
HNSC_metadata <- as.data.frame(colData(data_HNSC))
saveRDS(HNSC_metadata,
          file = "data/other_tumors/HNSC/TCGA-HNSC_546_samples_metadata.RDS")

##### OVARIAN CANCER ######

query_OV <- GDCquery(project = "TCGA-OV",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_OV,
            method = "api",
            files.per.chunk = 10)
data_OV <- GDCprepare(query_OV)

# Raw count data.
OV_count_data <- as.data.frame(assay(data_OV))
saveRDS(OV_count_data,
          file = "data/other_tumors/OV/TCGA-OV_379_samples_raw_count_data.RDS")

# Metadata.
OV_metadata <- as.data.frame(colData(data_OV))
saveRDS(OV_metadata,
          file = "data/other_tumors/OV/TCGA-OV_379_samples_metadata.RDS")

##### PANCREATIC CANCER ######

query_PAAD <- GDCquery(project = "TCGA-PAAD",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts",
                     experimental.strategy = "RNA-Seq",
                     legacy = FALSE)
GDCdownload(query_PAAD,
            method = "api",
            files.per.chunk = 10)
data_PAAD <- GDCprepare(query_PAAD)

# Raw count data.
PAAD_count_data <- as.data.frame(assay(data_PAAD))
saveRDS(PAAD_count_data,
          file = "data/other_tumors/PAAD/TCGA-PAAD_182_samples_raw_count_data.RDS")

# Metadata.
PAAD_metadata <- as.data.frame(colData(data_PAAD))
saveRDS(PAAD_metadata,
          file = "data/other_tumors/PAAD/TCGA-PAAD_182_samples_metadata.RDS")


##### UTERINE CANCER ######

query_UCEC <- GDCquery(project = "TCGA-UCEC",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts",
                     experimental.strategy = "RNA-Seq",
                     legacy = FALSE)
GDCdownload(query_UCEC,
            method = "api",
            files.per.chunk = 10)
data_UCEC <- GDCprepare(query_UCEC)

# Raw count data.
UCEC_count_data <- as.data.frame(assay(data_UCEC))
saveRDS(UCEC_count_data,
          file = "data/other_tumors/UCEC/TCGA-UCEC_587_samples_raw_count_data.RDS")

# Metadata.
UCEC_metadata <- as.data.frame(colData(data_UCEC))
saveRDS(UCEC_metadata,
          file = "data/other_tumors/UCEC/TCGA-UCEC_587_samples_metadata.RDS")

##### BLADDER CANCER ######

query_BLCA <- GDCquery(project = "TCGA-BLCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_BLCA,
            method = "api",
            files.per.chunk = 10)
data_BLCA <- GDCprepare(query_BLCA)

# Raw count data.
BLCA_count_data <- as.data.frame(assay(data_BLCA))
saveRDS(BLCA_count_data,
          file = "data/other_tumors/BLCA/TCGA-BLCA_433_samples_raw_count_data.RDS")

# Metadata.
BLCA_metadata <- as.data.frame(colData(data_BLCA))
saveRDS(BLCA_metadata,
          file = "data/other_tumors/BLCA/TCGA-BLCA_433_samples_metadata.RDS")

##### GALLBLADDER CANCER ######

query_CHOL <- GDCquery(project = "TCGA-CHOL",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_CHOL,
            method = "api",
            files.per.chunk = 10)
data_CHOL <- GDCprepare(query_CHOL)

# Raw count data.
CHOL_count_data <- as.data.frame(assay(data_CHOL))
saveRDS(CHOL_count_data,
          file = "data/other_tumors/CHOL/TCGA-CHOL_45_samples_raw_count_data.RDS")

# Metadata.
CHOL_metadata <- as.data.frame(colData(data_CHOL))
saveRDS(CHOL_metadata,
          file = "data/other_tumors/CHOL/TCGA-CHOL_45_samples_metadata.RDS")

##### ADRENAL GLAND CANCER ######

query_ACC <- GDCquery(project = "TCGA-ACC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_ACC,
            method = "api",
            files.per.chunk = 10)
data_ACC <- GDCprepare(query_ACC)

# Raw count data.
ACC_count_data <- as.data.frame(assay(data_ACC))
saveRDS(ACC_count_data,
          file = "data/other_tumors/ACC/TCGA-ACC_79_samples_raw_count_data.RDS")

# Metadata.
ACC_metadata <- as.data.frame(colData(data_ACC))
saveRDS(ACC_metadata,
          file = "data/other_tumors/ACC/TCGA-ACC_79_samples_metadata.RDS")

##### PHEOCHROMOCYTOMA AND PARAGANGLIOMA ######

query_PCPG <- GDCquery(project = "TCGA-PCPG",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      experimental.strategy = "RNA-Seq",
                      legacy = FALSE)
GDCdownload(query_PCPG,
            method = "api",
            files.per.chunk = 10)
data_PCPG <- GDCprepare(query_PCPG)

# Raw count data.
PCPG_count_data <- as.data.frame(assay(data_PCPG))
saveRDS(PCPG_count_data,
          file = "data/other_tumors/PCPG/TCGA-PCPG_186_samples_raw_count_data.RDS")

# Metadata.
PCPG_metadata <- as.data.frame(colData(data_PCPG))
saveRDS(PCPG_metadata,
          file = "data/other_tumors/PCPG/TCGA-PCPG_186_samples_metadata.RDS")

##### MELANOMA ######

query_SKCM <- GDCquery(project = "TCGA-SKCM",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts",
                       experimental.strategy = "RNA-Seq",
                       legacy = FALSE)
GDCdownload(query_SKCM,
            method = "api",
            files.per.chunk = 10)
data_SKCM <- GDCprepare(query_SKCM)

# Raw count data.
SKCM_count_data <- as.data.frame(assay(data_SKCM))
saveRDS(SKCM_count_data,
          file = "data/other_tumors/SKCM/TCGA-SKCM_472_samples_raw_count_data.RDS")

# Metadata.
SKCM_metadata <- as.data.frame(colData(data_SKCM))
saveRDS(SKCM_metadata,
          file = "data/other_tumors/SKCM/TCGA-SKCM_472_samples_metadata.RDS")