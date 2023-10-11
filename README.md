# Glioblastoma Discovery

<p align="center">
  <img src="https://github.com/marilenehohmuth/GlioblastomaDiscovery/blob/main/images/logo.png?raw=true" height="200"/>
</p>

## Overview

This repo provides the scripts used for the "Transcriptomics reveals an enhanced association between prion protein gene expression and vesicle dynamics signatures in glioblastomas" paper. The scripts are organized in the scripts folder and perform the paper analysis in the following manner: 

* single-cell_analyses.R: analysis of the Darmanis et al, Neftel et al, and Richards et al single-cell datasets.
* single-cell_protein_classes.R: correlations with the Panther protein classes.
* TCGA-GBM_analysis.R: analysis of GBM data from TCGA.
* TCGA-GBM_43_genes_analysis.R: analysis of GBM data from TCGA restricted to selected genes.
* TCGA_survival_analyses.R: survival analysis of GBM data from TCGA.
* TCGA_other_tumors_download.R: download of TCGA data for other tumors.
* TCGA_other_tumors_analyses.R: analysis of TCGA data for other tumors.

The repo also includes the code for the Shiny App in the shiny_app folder.
