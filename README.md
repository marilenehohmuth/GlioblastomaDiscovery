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

## Warning

Due to storage limitations we couldn't provide all genes to our database. If you try to run a gene which isn't in our database the app will not analyse it. With this is mind, the best way to check our avaliable genes is seeking for genes in the correlation tables.

## Tutorial
### Installing git 

Firstly, install git in your machine according with its Operational System:

#### Linux

Press CTRL + T in order to open your terminal, subsequently copy and paste the following code:

If you're on a Fedora distribution,

```
sudo dnf install git-all
```

If you're on a Debian distribution,

```
sudo apt install git-all
```

#### MacOs

Press Command + Space and type in "Terminal" and then open it. On the terminal copy and paste the following code:

If you have Homebrew,

```
brew install git
```

If you have MacPorts,

```
brew install git
```

(In case of not having any of them, check for documentation at [Homebrew](https://docs.brew.sh/) and [MacPorts](https://guide.macports.org/))

#### Windows

Use the installer available at [Git for Windows](https://git-scm.com/download/win).

### Cloning this repository

Now that you have Git installed, open your terminal and copy and paste the following code:


```
git clone https://github.com/marilenehohmuth/GlioblastomaDiscovery
```

This will make a copy of all files of this project in your computer.

### Installing R and RStudio

In order to run the code and have acess to an user-friendly interface you'll need to install R and RStudio, respectively.

#### R on Linux

Press CTRL + T in order to open your terminal, subsequently copy and paste the following code:

If you're on a Fedora distribution,

```
sudo dnf install R
```

If you're on a Debian distribution,

```
sudo apt install R
```

#### R on MacOS

Use the installer available at [R for MacOS](http://lib.stat.cmu.edu/R/CRAN/)

#### R on Windows

Use the installer available at [R for Windows](http://lib.stat.cmu.edu/R/CRAN/)

#### RStudio on Linux

At [RStudio Desktop](https://posit.co/download/rstudio-desktop/) search for the appropriate version based on your distribution and download it.

Press CTRL + T in order to open your terminal, subsequently copy and paste the following code:

If you're on a Fedora distribution,

```
sudo dnf install rstudio-desktop
```

If you're on a Debian distribution,

```
sudo apt install rstudio
```


#### RStudio on MacOS

Use the installer available at [RStudio for MacOS](https://posit.co/downloads/) 

#### RStudio on Windows

Use the installer available at [RStudio for Windows](https://posit.co/download/rstudio-desktop/)

After those steps you'll be able to visualize our app code at your own machine. However, in order to get access to the user-friendly interface and use the app tools, you need to install all required R libraries. Some of them are in the code and will be automatically installed as soon as you try to run the app. Some of them will need to be individually installed.
