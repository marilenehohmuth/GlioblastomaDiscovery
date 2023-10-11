# Glioblastoma

library(devtools)

# Install fusca dependencies

devtools::install_github("sctyner/geomnet")

list.of.packages <- c('cccd', 'grid', 'tsne', 'Rtsne', 'igraph', 'mclust', 'ggplot2', 'pheatmap', 'reshape', 'reshape2')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, repos = c("http://cran.rstudio.com/", "https://bioconductor.org/biocLite.R"))

install.packages("Vennerable", repos = "http://R-Forge.R-project.org")

# Install fusca
devtools::install_github("edroaldo/fusca")

# Load packages ---------------------------------------------------------------

library(shiny)
library(fusca)
library(dplyr)
library(ggplot2)
library(Matrix)

source('./functions/modified_functions.R')

# Process object ----------------------------------------------------------

#' Classify cells according to gene expression 
#'
#' Classify the cells expressing a gene into positive and negative cells.
#'
#' @param selected_gene character; name of the gene.
#' @param cellrouter CellRouter object; the object to be processed.
classify_cells <- function(cellrouter, selected_gene){
  # Get expression
  gene_exp <- cellrouter@assays$RNA@ndata[selected_gene,]
  gene_exp <- as.data.frame(gene_exp)
  colnames(gene_exp) <- 'GENE'
  # Classify
  gene_exp <- gene_exp %>% 
    mutate(classification = case_when(GENE == 0 ~ paste0(selected_gene, '_negative_cells'),
                                      GENE > 0 ~ paste0(selected_gene, '_positive_cells')))
  # Add classification.
  cellrouter <- addInfo(object = cellrouter,
                        metadata = gene_exp,
                        colname = paste0(selected_gene, '_status'),
                        metadata.column = "classification")
  return(cellrouter)
}

#' Calculate gene correlation 
#'
#' Perform gene correlation analysis.
#'
#' @param selected_gene character; name of the gene.
#' @param cellrouter CellRouter object; the object to be processed.
calculate_gene_correlation <- function(cellrouter, selected_gene){
  print("Starting correlation analysis.")
  df_darmanis <- data.frame(matrix(0, ncol=4))
  colnames(df_darmanis) <- c('gene1', 'gene2', 'correlation', 'pvalue')
  i=1
  # here ? correlation is calculated
  for(g in rownames(cellrouter@assays$RNA@ndata)[1:25]){
    c <- cor.test(as.numeric(cellrouter@assays$RNA@ndata[g,]), as.numeric(cellrouter@assays$RNA@ndata[selected_gene,]))
    df_darmanis[i,'gene1'] <- selected_gene
    df_darmanis[i,'gene2'] <- g
    df_darmanis[i,'correlation'] <- c$estimate #stores the correlation
    df_darmanis[i,'pvalue'] <- c$p.value #stores the correlation
    i <- i + 1
  }
  print(paste('Computed correlation between', selected_gene, 'and all genes.'))
  return(df_darmanis)
}

#' Process CellRouter 
#'
#' Applies the classify cells and calculate gene correlation functions to the
#' CellRouter object.
#'
#' @param selected_gene character; name of the gene.
#' @param cellrouter CellRouter object; the object to be processed.
process_cellrouter <- function(cellrouter, selected_gene){
  return(classify_cells(cellrouter, selected_gene))
}


# Datasets ---------------------------------------------------------------

path_darmanis <- './data/reduced_darmanis.RDS'
path_neftel <- './data/reduced_neftel.RDS'
path_richards <- './data/reduced_richards.RDS'

path_mutual_genes <- './data/overall_hvgs.RDS'

all_genes <- readRDS(path_mutual_genes)



# Read text and figure ----------------------------------------------------

intro_text_path <- "./design/intro_text.txt"
intro_text <- readChar(file.path(intro_text_path), file.info(intro_text_path)$size)
image_path <- file.path("./design/logo.png")
image_type <- "image/png"
darmanis_text_path <- "./design/darmanis_text.txt"
darmanis_text <- readChar(file.path(darmanis_text_path), file.info(darmanis_text_path)$size)
neftel_text_path <- "./design/neftel_text.txt"
neftel_text <- readChar(file.path(neftel_text_path), file.info(neftel_text_path)$size)
richards_text_path <- "./design/richards_text.txt"
richards_text <- readChar(file.path(richards_text_path), file.info(richards_text_path)$size)



# UI ----------------------------------------------------------------------

# Define UI
ui <- navbarPage("GBMdiscovery",
                 # Theme
                 # https://rstudio.github.io/shinythemes/
                 theme = bslib::bs_theme(bootswatch = 'lumen'),
                 # Feedback
                 # shinyFeedback::useShinyFeedback(),
                 
                 # Description -------------------------------------------------------------
                 tabPanel('Description',
                          fluidRow(
                            # Quick intro.
                            column(8, 
                                   h3('GBMdiscovery'),
                                   helpText(intro_text),
                                   imageOutput('photo')
                            ),
                            column(4, 
                                   textInput('input_gene', h3('Selected gene'), 
                                             value='PRNP'),
                                   actionButton('update_analysis', 
                                                'Update analysis', 
                                                class='btn btn-success',
                                                icon=icon('sync'))
                            )
                            
                          )
                 ),
                 
                 
                 
                 # Darmanis ------------------------------------------------------------------
                 
                 # Darmanis data
                 tabPanel('Darmanis',
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                            sidebarPanel(
                              h3('Darmanis'),
                              helpText(darmanis_text),
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              tabsetPanel(
                                tabPanel('Marker gene signatures',
                                         plotOutput('dar_rd')
                                ),
                                tabPanel('Gene expression',
                                         plotOutput('dar_ge')
                                ),
                                tabPanel('Correlations',
                                         dataTableOutput('dar_cor')
                                ) # ,
                                # tabPanel('Altered biological phenomena',
                                #          plotOutput('dar_rp')
                                # )
                              )
                            )
                          )
                 ),
                 
                 
                 # Neftel ------------------------------------------------------------------
                 
                 # - data
                 tabPanel('Neftel',
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                            sidebarPanel(
                              h3('Neftel'),
                              helpText(neftel_text),
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              tabsetPanel(
                                tabPanel('Marker gene signatures',
                                         plotOutput('nef_rd')
                                ),
                                tabPanel('Gene expression',
                                         plotOutput('nef_ge')
                                ),
                                tabPanel('Correlations',
                                         dataTableOutput('nef_cor')
                                ) # ,
                                # tabPanel('Altered biological phenomena',
                                #          plotOutput('nef_rp')
                                # )
                              )
                            )
                          )
                 ),
                 
                 
                 # Richards ----------------------------------------------------------------
                 
                 tabPanel('Richards',
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                            sidebarPanel(
                              h3('Richards'),
                              helpText(richards_text),
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              tabsetPanel(
                                tabPanel('Marker gene signatures',
                                         plotOutput('ric_rd')
                                ),
                                tabPanel('Gene expression',
                                         plotOutput('ric_ge')
                                ),
                                tabPanel('Correlations',
                                         dataTableOutput('ric_cor')
                                ) # ,
                                # tabPanel('Altered biological phenomena',
                                #          plotOutput('ric_rp')
                                # )
                              )
                            )
                          )
                 ),
                 ## removing data 4 (bulk RNA)
                 # tabPanel('Data 4',
                 #          # Sidebar with a slider input for number of bins 
                 # ),
                 ## removed the analysis for the combined datasets
                 # tabPanel('Combined Datasets',
                 #          # 
                 #          sidebarLayout(
                 #            sidebarPanel(
                 #              h3('Genes from combined datasets'),
                 #              helpText('Here we present the genes that are present in the combined datasets.'),
                 #            ),
                 #            # Show a plot of the generated distribution
                 #            mainPanel(
                 #              tabsetPanel(
                 #                tabPanel('Single-cell datasets',
                 #                         dataTableOutput('sc_genes')
                 #                ),
                 #                tabPanel('All datasets',
                 #                         dataTableOutput('all_genes')
                 #                )
                 #              )
                 #            )
                 #          )
                 # )
                 
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Reactive variables
  selected_gene <- reactiveVal('PRNP')
  cellrouter_darmanis <- reactiveVal(process_cellrouter(readRDS(path_darmanis), 'PRNP'))
  cellrouter_neftel <- reactiveVal(process_cellrouter(readRDS(path_neftel), 'PRNP'))
  cellrouter_richards <- reactiveVal(process_cellrouter(readRDS(path_richards), 'PRNP'))
  correlation_darmanis <- reactiveVal(calculate_gene_correlation(readRDS(path_darmanis), 'PRNP'))
  correlation_neftel <- reactiveVal(calculate_gene_correlation(readRDS(path_neftel), 'PRNP'))
  correlation_richards <- reactiveVal(calculate_gene_correlation(readRDS(path_richards), 'PRNP'))
  
  # Gene name
  observeEvent(input$update_analysis, {
    # Message is not presented.
    validate(
      need(input$input_gene %in% all_genes, 'Provide gene in list')
    )
    selected_gene(input$input_gene)
    cellrouter_darmanis(process_cellrouter(cellrouter_darmanis(), selected_gene()))
    cellrouter_neftel(process_cellrouter(cellrouter_neftel(), selected_gene()))
    cellrouter_richards(process_cellrouter(cellrouter_richards(), selected_gene()))
    correlation_darmanis(calculate_gene_correlation(cellrouter_darmanis(), selected_gene()))
    correlation_neftel(calculate_gene_correlation(cellrouter_neftel(), selected_gene()))
    correlation_richards(calculate_gene_correlation(cellrouter_richards(), selected_gene()))
  })
  
  # Description
  output$photo <- renderImage({
    list(
      src = image_path,
      contentType = image_type,
      width = '60%'
    )
  }, deleteFile = FALSE)
  
  # Darmanis ----------------------------------------------------------------
  # Darmanis
  output$dar_rd <- renderPlot({plotReducedDimension(cellrouter_darmanis(),
                                                    reduction.type='umap',
                                                    annotation=paste0(selected_gene(), '_status'), 
                                                    annotation.color=paste0(selected_gene(), '_status_color'),
                                                    showlabels=FALSE,
                                                    dotsize=1, 
                                                    labelsize=5, 
                                                    convex=FALSE)
  })
  output$dar_ge <- renderPlot({plotDRExpression(cellrouter_darmanis(),
                                                assay.type = 'RNA',
                                                genelist = c(selected_gene()),
                                                reduction.type = 'umap',
                                                threshold = 2,
                                                dotsize = 1,
                                                title = "")})
  output$dar_cor <- renderDataTable(calculate_gene_correlation(cellrouter_darmanis(), selected_gene()), 
                                    options = list(pageLength = 5))
  
  # Neftel ------------------------------------------------------------------
  output$nef_rd <- renderPlot({plotReducedDimension(cellrouter_neftel(),
                                                    reduction.type='umap',
                                                    annotation=paste0(selected_gene(), '_status'), 
                                                    annotation.color=paste0(selected_gene(), '_status_color'),
                                                    showlabels=FALSE,
                                                    dotsize=1, 
                                                    labelsize=5, 
                                                    convex=FALSE)
  })
  output$nef_ge <- renderPlot({plotDRExpression(cellrouter_neftel(),
                                                assay.type = 'RNA',
                                                genelist = c(selected_gene()),
                                                reduction.type = 'umap',
                                                threshold = 2,
                                                dotsize = 1,
                                                title = "")})
  output$nef_cor <- renderDataTable(calculate_gene_correlation(cellrouter_neftel(), selected_gene()), 
                                    options = list(pageLength = 5))
  
  # Richards ----------------------------------------------------------------
  output$ric_rd <- renderPlot({plotReducedDimension(cellrouter_richards(),
                                                    reduction.type='umap',
                                                    annotation=paste0(selected_gene(), '_status'), 
                                                    annotation.color=paste0(selected_gene(), '_status_color'),
                                                    showlabels=FALSE,
                                                    dotsize=1, 
                                                    labelsize=5, 
                                                    convex=FALSE)
  })
  output$ric_ge <- renderPlot({plotDRExpression(cellrouter_richards(),
                                                assay.type = 'RNA',
                                                genelist = c(selected_gene()),
                                                reduction.type = 'umap',
                                                threshold = 2,
                                                dotsize = 1,
                                                title = "")})
  output$ric_cor <- renderDataTable(calculate_gene_correlation(cellrouter_richards(), selected_gene()), 
                                    options = list(pageLength = 5))
  
  
  
  # Combined ----------------------------------------------------------------
  
  # Combined
  output$sc_genes <- renderDataTable(data.frame(Darmanis = rownames(cellrouter_darmanis()@assays$RNA@ndata),
                                                Neftel = rownames(cellrouter_neftel()@assays$RNA@ndata),
                                                Richards = rownames(cellrouter_richards()@assays$RNA@ndata)), 
                                     options = list(pageLength = 5,
                                                    dom = 'ft')
  ) 
  output$all_genes <- renderDataTable(data.frame('All genes' = all_genes), 
                                      options = list(pageLength = 5))
  
}

# Run the application 
shinyApp(ui = ui, server = server)