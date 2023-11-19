
# ANALYSIS OF SINGLE-CELL DATASETS ----------------------------------------

# @ This script generates the plots present in the second main figure of the # 
# @ manuscript + supplemental material.                                      #

## Loading packages --------------------------------------------------------

library(fusca)                    # 1.2.1
library(ggplot2)                  # 3.3.5
library(dplyr)                    # 1.0.7
library(Matrix)                   # 1.3-4
library(RColorBrewer)             # 1.1-2
library(Hmisc)                    # 4.5-0
library(corrplot)                 # 0.92
library(stringr)                  # 1.5.0
library(doParallel)               # 1.0.17
library(clusterProfiler)          # 4.8.1
library(org.Hs.eg.db)             # 3.17.0
library(irlba)                    # 2.3.5.1
library(VennDiagram)              # 1.7.3

## Defining some functions -------------------------------------------------

# Load sample information table.
sample_info <- read.csv(
    paste0(getwd(), "/data/sample_info.csv"),
    header = TRUE,
    row.names = NULL
)
filter_out_samples <- function(
    mdata, # A metadata table where rows are cells and columns are cell metadata.
    sample_column # Name of column in metadata table that contains sample IDs.
) {

    # Keep only samples marked as "YES" in sample information table.
    mdata <- mdata[mdata[,sample_column] %in% sample_info$sample[sample_info$use_in_analyses == "YES"],]

    # Return filtered metadata table.
    return(mdata)
}

# process_cellrouter()
# @ This function processes a given CellRouter object & saves the processed object to an output file.
process_cellrouter <- function(
    cellrouter, # CellRouter object.
    metadata, # Metadata table.
    normalized, # Boolean indicating if count data is already normalized.
    sample_column_id, # Column that contains sample IDs in metadata table.
    dataset # Dataset name.
) {

    cellrouter@assays$RNA@rawdata <- cellrouter@assays$RNA@rawdata[rownames(cellrouter@assays$RNA@ndata), colnames(cellrouter@assays$RNA@ndata)]

    # Add sample IDs to the CellRouter object.
    cellrouter@assays$RNA@sampTab$samples <- metadata[,sample_column_id][match(rownames(cellrouter@assays$RNA@sampTab), rownames(metadata))]

    # Normalize count data if not yet normalized.
    if(normalized == FALSE) {
        cellrouter <- Normalize(cellrouter)
    } else {
        cellrouter@assays$RNA@ndata <- cellrouter@assays$RNA@rawdata
    }

    # Scale data according to all genes.
    cellrouter <- scaleData(cellrouter, blocksize = nrow(cellrouter@assays$RNA@ndata))

    # Perform PCA.
    cellrouter <- computePCA(cellrouter, num.pcs = 50, seed = 42)

    # Compute UMAP.
    cellrouter <- computeUMAP(cellrouter, num.pcs = 15, seed = 42)

    # Save processed CellRouter object.
    save(cellrouter, file = paste0("results/", dataset, "/", str_to_title(dataset), "_Normalized_Scaled_PCA_UMAP_CellRouter_object.R"))

    # Return processed CellRouter object.
    return(cellrouter)
}

# plot_prnp_across_samples()
# @ This function plots on a boxplot the PRNP expression across samples in a given dataset.
plot_prnp_across_samples <- function(
    cellrouter, # CellRouter object.
    dataset, # Dataset name.
    plot_width, # Boxplot width.
    plot_height # Boxplot height
) {

    # Transpose normalized data & ensure row names have the same order as row names in the sampTab slot of the RNA assay.
    transposed_norm_data <- as.data.frame(t(as.matrix(cellrouter@assays$RNA@ndata)))
    transposed_norm_data <- transposed_norm_data[rownames(cellrouter@assays$RNA@sampTab),]

    # Concatenate dataframes to have gene expression and metadata information in the same table.
    expr <- cbind(transposed_norm_data, cellrouter@assays$RNA@sampTab)

    # Create boxplot.
    plot <- ggplot(
        expr,
        aes(x = reorder(samples, PRNP), y = PRNP)
    ) +
        geom_violin(aes(fill = samples), scale = "width") +
        geom_boxplot(fill = "white", outlier.shape = NA, width = 0.1) +
        xlab("Sample") +
        ylab(expression("Normalized"~italic("PRNP")~"expression")) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 12)
        )
    
    # Save plot to output file.
    pdf(
        paste0("results/", dataset, "/", str_to_title(dataset), "_PRNP_expression_across_samples_boxplot.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()
    
}

# add_prnp_status()
# @ This function takes a CellRouter object and assigns each cell a PRNP status.
add_prnp_status <- function(cellrouter) {
    
    # Getting cells that express PRNP.
    prnp_pos <- colnames(cellrouter@assays$RNA@ndata)[cellrouter@assays$RNA@ndata["PRNP",] > 0]
    # Getting cells that do not express PRNP.
    prnp_neg <- colnames(cellrouter@assays$RNA@ndata)[cellrouter@assays$RNA@ndata["PRNP",] == 0]

    # Adding PRNP status.
    cellrouter@assays$RNA@sampTab <- cellrouter@assays$RNA@sampTab %>% mutate(PRNP_status = case_when(
        rownames(cellrouter@assays$RNA@sampTab) %in% prnp_pos ~ "PRNP_positive_cells",
        rownames(cellrouter@assays$RNA@sampTab) %in% prnp_neg ~ "PRNP_negative_cells")
    )

    # Return processed CellRouter object.
    return(cellrouter)
}

# concatenate_data()
# @ This function puts together cell embeddings, gene expression and cell metadata.
concatenate_data <- function(cellrouter) {

    # Transpose normalized data to get cells as row names.
    transposed_norm_data <- as.data.frame(t(as.matrix(cellrouter@assays$RNA@ndata)))

    # Ensure row names in cell embeddings of UMAP slot have the same order as the row names of the transposed, normalized count data.
    cellrouter@umap$cell.embeddings <- cellrouter@umap$cell.embeddings[rownames(transposed_norm_data),]

    # Concatenate dataframes to have cell embeddings and gene expression information in the same table.
    dim.red <- cbind(cellrouter@umap$cell.embeddings, transposed_norm_data)
    colnames(dim.red)[1] <- "UMAP1"
    colnames(dim.red)[2] <- "UMAP2"

    # Ensure row names have the same order as the row names in the sampTab slot of the RNA assay.
    dim.red <- dim.red[rownames(cellrouter@assays$RNA@sampTab),]

    # Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
    dim.red <- cbind(dim.red, cellrouter@assays$RNA@sampTab)

    # Return dataframe.
    return(dim.red)
}

# plot_prnp_umap()
# @ This function projects cells onto a UMAP plot and color them according to PRNP expression, saving the plot to an output file.
plot_prnp_umap <- function(
    dim.red, # Dataframe with cell embeddings, gene expression and cell metadata.
    dataset, # Dataset name.
    plot_width, # Width for UMAP plot.
    plot_height # Height for UMAP plot.
) {

    plot <- ggplot(
        dim.red %>% dplyr::arrange(PRNP),
        aes(x = UMAP1, y = UMAP2, color = PRNP)
    ) +
        geom_point(size = 0.1) +
        theme(
            legend.position = "bottom",
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
            legend.key.width = unit(1.2, 'cm'),
            legend.key.height = unit(0.5, 'cm')
        ) +
        scale_color_gradientn(
            colors = colorRampPalette((brewer.pal(9, "RdPu")))(1000),
            limits = c(min(dim.red$PRNP), max(dim.red$PRNP)),
            name = expression(~italic("PRNP")~"expression  ")
        )

    # Save plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_PRNP_expression_across_cells_UMAP_plot.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()

}

# plot_sample_umap()
# @ This function projects cells onto a UMAP plot and color them according to sample IDs, saving the plot to an output file.
plot_sample_umap <- function(
    dim.red, # Dataframe with cell embeddings, gene expression and cell metadata.
    dataset, # Dataset name.
    plot_width, # Width for UMAP plot.
    plot_height # Height for UMAP plot.
) {
    # Get number of samples.
    n_samples <- dim.red$samples %>% unique %>% length

    # Define number of columns for plot legend.
    n_cols <- 2 
    if(n_samples >= 15) n_cols <- 8

    plot <- ggplot(
        dim.red,
        aes(x = UMAP1, y = UMAP2, color = samples)
    ) +
        geom_point(size = 0.1) +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size=5),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
            legend.key.size = unit(2, 'cm')
        ) + 
        guides(color = guide_legend(
            override.aes = list(size=0.8),
            keywidth = 0,
            keyheight = 0,
            ncol = n_cols,
            title = "Sample")
        )

    # Save plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_samples_UMAP_plot.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()

}

# get_signatures()
# @ This function uses CellRouter's findSignatures() function to find the marker genes of PRNP+ and PRNP- subpopulations.
get_signatures <- function(
    cellrouter, # CellRouter object.
    dataset # Dataset name.
) {

    # Find PRNP+ and PRNP- marker gene signatures & save results to output file.
    markers <- findSignatures(
        object = cellrouter,
        assay.type = "RNA",
        column = "PRNP_status",
        test.use = "wilcox",
        pos.only = FALSE,
        fc.threshold = 0
    )

    # Save results to output file.
    write.csv(
        markers, 
        file = paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_PRNP+_vs_PRNP-_signature_marker_genes.csv")
    )

    # Return marker genes.
    return(markers)
}

# plot_signature_volcano()
# @ This function uses get_signatures()'s output to create a volcano plot.
plot_signature_volcano <- function(
    markers, # Dataframe with marker genes as output by CellRouter's findSignatures().
    plot_width, # Width for volcano plot.
    plot_height # Height for volcano plot.
) {

    # Create volcano plot.
    plot <- ggplot(
        markers,
        aes(x = fc, y = -log10(pval))
    ) +
    geom_point(size = 1, fill = "royalblue", shape = 21) +
    theme_classic() + 
    theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
        ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ylab(expression(-log[10](p-value))) +
    xlab(expression(log[2]("Fold change")))

    # Save plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_Correlation_PRNPvsAllGenes_plot.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()

}

# get_correlations_with_prnp()
# @ This function calculates the Pearson correlation between PRNP expression and the expression of each gene based on normalized counts
# of a given dataset.
get_correlations_with_prnp <- function(
    cellrouter, # CellRouter object.
    dataset # Dataset name.
) {

    # Create dataframe to store correlations.
    df <- data.frame(matrix(0, ncol = 4))
    colnames(df) <- c('gene1', 'gene2', 'correlation', 'pvalue')

    # Iterate over all genes present in the ndata slot of the RNA assay, calculating the Pearson correlation between their
    # expression & PRNP expression.
    i <- 1
    for(gene in rownames(cellrouter@assays$RNA@ndata)){
        c <- cor.test(
            as.numeric(cellrouter@assays$RNA@ndata[gene,]), 
            as.numeric(cellrouter@assays$RNA@ndata['PRNP',]),
            method = "pearson")
        df[i,'gene1'] <- 'PRNP'
        df[i,'gene2'] <- gene
        df[i,'correlation'] <- c$estimate # Store the correlation.
        df[i,'pvalue'] <- c$p.value # Store the p-value.
        i <- i + 1
    }
    # Save complete results to output file.
    write.csv(df, file = paste0("results/", dataset, "/", str_to_title(dataset), "_Correlation_PRNPvsAllGenes.csv"))

    # Remove cases for which correlation is NA.
    df <- df[!is.na(df$correlation),]

    # Reorder the dataframe according to correlation coeffiecients in decreasing order.
    df <- df[order(df$correlation, decreasing = TRUE),]

    # Remove correlation PRNP vs PRNP (which is equal to 1).
    df <- df[df$gene2 != "PRNP",]

    # Add row names to dataframe that will determine the position of the genes later on in the correlation plot.
    rownames(df) <- 1:length(rownames(df))

    # Add classifications.
    df <- df %>% mutate(status = case_when(
        correlation < 0 & pvalue <= 0.05 ~ "Negative correlation",
        correlation > 0 & pvalue <= 0.05 ~ "Positive correlation",
        correlation == 0 ~ "Non-significant",
        pvalue > 0.05 ~ "Non-significant"
    ))

    # Save filtered & classified results to output file.
    write.csv(
        df, 
        file = paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_Correlation_PRNPvsAllGenes_filt_withStatus.csv")
    )

    # Return dataframe with correlations.
    return(df)
}

# plot_correlations_with_prnp()
# @ This function uses get_correlations_with_prnp()'s output to create a correlation plot.
plot_correlations_with_prnp <- function(
    df, # Dataframe as output by get_correlations_with_prnp().
    dataset, # Dataset name.
    plot_width, # Width for correlation plot.
    plot_height # Height for correlation plot.
) {

    # Create plot.
    plot <- ggplot(
        df,
        aes(x = as.numeric(rownames(df)), y = as.numeric(correlation), color = status)
    ) +
        geom_point() +
        theme_bw() +
        geom_hline(yintercept = 0, linetype = 'dotted') +
        xlab("Ranked genes") +
        ylab("Pearson correlation") +
        scale_color_manual(values = c("royalblue1", "snow3", "indianred1"), name = "Status")  +
        ggtitle(paste0(str_to_title(dataset), " et al")) +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)
        )

    # Save plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_Correlation_PRNPvsAllGenes_plot.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()

}

# do_gsea()
# @ This function uses findSignatures()'s output to perform GSEA & save results to output csv file.
do_gsea <- function(
    df, # Dataframe containing findSignatures results.
    dataset, # Dataset name.
    population_ # Cell population name, as output by findSignatures.
) {
    # Select only cell population specified.
    df <- df[df$population == population_,]

    # Filter out non-significant results.
    df <- df[df$pval <= 0.05,]

    # Get ranked gene list.
    gene_list <- df$fc
    names(gene_list) <- df$gene
    gene_list <- gene_list[!is.na(gene_list)]
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Perform Gene Set Enrichment Analysis (GSEA) in the specified population.
    gse <- gseGO(
        geneList = gene_list, 
        ont = "ALL", 
        keyType = "SYMBOL", 
        nPerm = 10000, 
        minGSSize = 3, 
        maxGSSize = 800, 
        pvalueCutoff = 0.05, 
        verbose = TRUE, 
        OrgDb = org.Hs.eg.db, 
        pAdjustMethod = "BH"
    )   
    gse <- as.data.frame(gse@result)

    # Save GSEA results table to output file.
    write.csv(
        gse, 
        file = paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_", population_, "_GSEA_gseGO_ontALL_padj0.05.csv")
    )

    # Return GSEA results table.
    return(gse)
}

# plot_gsea()
# @ This function uses gseGO()'s output to plot GSEA results, highlighting selected terms.
plot_gsea <- function(
    gse, # Dataframe containing gseGO results.
    dataset, # Dataset name.
    population_, # Cell population name used to perform GSEA with gseGO.
    title, # Title for plot.
    terms, # Terms to highlight on the plot.
    plot_width, # Width for correlation plot.
    plot_height # Height for correlation plot.
) {
    # Create dot plot with GSEA results.
    plot <- ggplot(
        gse[gse$Description %in% terms,],
        aes(x = NES, y = reorder(Description, NES), fill = -log10(p.adjust), size = setSize) 
    ) +
        geom_point(shape = 21, color = "black") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        ) +
        xlab("Normalized enrichment scores (NES)") +
        ylab("Gene Ontology (GO) term") +
        ggtitle(title) +
        scale_colour_gradientn(
            colors = colorRampPalette(rev(brewer.pal(9, "Reds")))(100),
            name = expression(-log[10]~"(Adjusted p-value)")
        ) +
        scale_size_continuous(range = c(1,5), name = "# Genes")

    # Save dot plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_", population_, "_GSEA_gseGO_ontALL_padj0.05_selectedTerms.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()
}

# do_ora()
# @ This function uses findSignatures()'s output to perform ORA & save results to output csv files.
do_ora <- function(
    df, # Dataframe containing findSignatures results.
    dataset, # Dataset name.
    population_ # Cell population name, as output by findSignatures.
) {
    # Select only cell population specified.
    df <- df[df$population == population_,]

    # Filter out non-significant results.
    df <- df[df$pval <= 0.05,]

    # Get gene lists.
    gene_list <- list(
        "Upregulated" = df$gene[df$fc > 0],
        "Downregulated" = df$gene[df$fc < 0],
    )

    # Perform Over-Representation Analysis (ORA) in the specified population.
    comparison <- compareCluster(
        gene_list,
        fun = "enrichGO",
        ont = "ALL",
        keyType = "SYMBOL",
        OrgDb = "org.Hs.eg.db",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
    )
    comparison <- as.data.frame(comparison)

    # Save ORA results table to output file.
    write.csv(
        comparison, 
        file = paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_", population_, "_ORA_compareCluster_enrichGO_ontALL_padj0.05.csv")
    )

    # Return ORA results.
    return(comparison)
}
    
# plot_ora()
# @ This function uses enrichGO()'s output to plot ORA results, highlighting selected terms for the selected population.
plot_ora <- function(
    comparison, # Dataframe containing enrichGO results generated with compareCluster.
    dataset, # Dataset name.
    population_, # Cell population name used to perform ORA with enrichGO.
    title, # Title for plot.
    terms, # Terms to highlight on the plot.
    plot_width, # Width for correlation plot.
    plot_height # Height for correlation plot.
) {
    # Create dot plot with ORA results.
    plot <- ggplot(
          comparison[comparison$Description %in% terms,],
          aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), fill = Cluster, size = Count) 
    ) +
        geom_point(shape = 21, color = "black") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        ) +
        xlab(expression(-log[10]~"(Adjusted p-value)")) +
        ylab("Gene Ontology (GO) term") +
        ggtitle(title) +
        scale_fill_manual(
            values = c("indianred1", "royalblue1"),
            labels = c("Upregulated", "Downregulated"),
            name = "Cluster"
        ) +
        scale_size_continuous(range = c(1,5), name = "# Genes") +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed")

    # Save dot plot to output file.
    pdf(
        paste0(getwd(), "/results/", dataset, "/", str_to_title(dataset), "_", population_, "_ORA_compareCluster_enrichGO_ontALL_padj0.05_selectedTerms.pdf"),
        width = plot_width,
        height = plot_height
    )
    print(plot)
    dev.off()
}


# Darmanis et al. dataset -------------------------------------------------


# Data downloaded from:
# http://gbmseq.org


## Step 1: Loading and subsetting data for analysis ------------------------

# Load Darmanis et al count data.
data_darmanis <- read.csv(
    paste0(getwd(), "/data/darmanis/GBM_raw_gene_counts.csv"),
    header = TRUE,
    row.names = 1,
    sep = " ",
    check.names = FALSE
)

# Load Darmanis et al metadata.
metadata_darmanis <- read.csv(
    paste0(getwd(), "/data/darmanis/GBM_metadata.csv"),
    header = TRUE,
    row.names = 1,
    sep = " "
)

# Subset metadata to keep only malignant cells.
malignant_metadata_darmanis <- metadata_darmanis[metadata_darmanis$Cluster_2d %in% c(1, 4, 11),]

# Filter out samples that should not be used for downstream analyses.
malignant_metadata_darmanis <- filter_out_samples(
    mdata = malignant_metadata_darmanis, 
    sample_column = "Sample.name"
)

# Subset count data to keep only malignant cells & cells from samples that should be kept for
# downstream analyses.
malignant_data_darmanis <- data_darmanis %>% dplyr::select(rownames(malignant_metadata_darmanis))

## Step 2: Data processing -------------------------------------------------

# Create CellRouter object for Darmanis et al dataset.
cellrouter_darmanis <- CreateCellRouter(
    as.matrix(malignant_data_darmanis),
    min.genes = 0,
    min.cells = 0,
    is.expr = 0
)

# Process Darmanis' CellRouter object.
cellrouter_darmanis <- process_cellrouter(
    cellrouter = cellrouter_darmanis,
    metadata = malignant_metadata_darmanis,
    normalized = FALSE,
    sample_column_id = "Sample.name",
    dataset = "darmanis"
)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
plot_prnp_across_samples(
    cellrouter = cellrouter_darmanis,
    dataset = "darmanis",
    plot_width = 2.5,
    plot_height = 3
)

## Step 3: Establishing groups with distinct PRNP levels -------------------

# Classifying cells as PRNP+ or PRNP-.
cellrouter_darmanis <- add_prnp_status(cellrouter_darmanis)

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.darmanis <- concatenate_data(cellrouter = cellrouter_darmanis)

#### @ FIGURE 2B, LEFT PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
plot_prnp_umap(
    dim.red = dim.red.darmanis,
    dataset = "darmanis",
    plot_width = 4.5,
    plot_height = 4
)

#### @ FIGURE S2A, LEFT PANEL (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample.
plot_sample_umap(
    dim.red = dim.red.darmanis,
    dataset = "darmanis",
    plot_width = 4.5,
    plot_height = 4
)

## Step 4: Finding signatures of PRNP+ and PRNP- cells ---------------------

markers_darmanis <- get_signatures(
    cellrouter = cellrouter_darmanis, 
    dataset = "darmanis"
)

# Currently not used in the manuscript.
# plot_signature_volcano(markers_darmanis, "darmanis")

# Perform GSEA on markers found for PRNP positive cells.
gsea_darmanis_prnp_pos <- do_gsea(
    df = markers_darmanis,
    dataset = "darmanis",
    population_ = "PRNP_positive_cells"
)
plot_gsea(
    gse = gsea_darmanis_prnp_pos,
    dataset = "darmanis",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform GSEA on markers found for PRNP negative cells.
gsea_darmanis_prnp_neg <- do_gsea(
    df = markers_darmanis,
    dataset = "darmanis",
    population_ = "PRNP_negative_cells"
)
plot_gsea(
    gse = gsea_darmanis_prnp_neg,
    dataset = "darmanis",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP positive cells.
ora_darmanis_prnp_pos <- do_ora(
    df = markers_darmanis,
    dataset = "darmanis",
    population_ = "PRNP_positive_cells"
)
plot_ora(
    gse = ora_darmanis_prnp_pos,
    dataset = "darmanis",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP negative cells.
ora_darmanis_prnp_neg <- do_ora(
    df = markers_darmanis,
    dataset = "darmanis",
    population_ = "PRNP_negative_cells"
)
plot_ora(
    gse = ora_darmanis_prnp_neg,
    dataset = "darmanis",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

## Step 5: Correlation analysis --------------------------------------------

df_darmanis <- get_correlations_with_prnp(
    cellrouter = cellrouter_darmanis, 
    dataset = "darmanis"
)

#### @ FIGURE 2C, LEFT PANEL (MAIN) @ #### 
plot_correlations_with_prnp(
    df = df_darmanis, 
    dataset = "darmanis",
    plot_width = 5,
    plot_height = 3
)


# Neftel et al. dataset ---------------------------------------------------


# Data downloaded from:
# https://singlecell.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma


## Step 1: Loading and subsetting data for analysis ------------------------

# Load Neftel et al count data.
data_neftel <- read.table(
    paste0(getwd(), "/data/neftel/IDHwtGBM.processed.SS2.logTPM.txt.gz"),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

# Load Neftel et al metadata.
metadata_neftel <- read.table(
    paste0(getwd(), "/data/neftel/IDHwt.GBM.Metadata.SS2.txt"),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

# Subset metadata to keep only malignant cells.
malignant_metadata_neftel <- metadata_neftel[metadata_neftel$CellAssignment == "Malignant",]

# Filter out samples that should not be used for downstream analyses.
malignant_metadata_neftel <- filter_out_samples(
    mdata = malignant_metadata_neftel, 
    sample_column = "Sample"
)

# Subset count data to keep only malignant cells & cells from samples that should be kept for
# downstream analyses.
malignant_data_neftel <- data_neftel %>% dplyr::select(rownames(malignant_metadata_neftel))

## Step 2: Data processing -------------------------------------------------

# Create CellRouter object for Neftel et al dataset.
cellrouter_neftel <- CreateCellRouter(
    as.matrix(malignant_data_neftel),
    min.cells = 0,
    min.genes = 0,
    is.expr = 0
)

# Process Neftel's CellRouter object.
cellrouter_neftel <- process_cellrouter(
    cellrouter = cellrouter_neftel,
    metadata = malignant_metadata_neftel,
    normalized = TRUE,
    sample_column_id = "Sample",
    dataset = "neftel"
)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
plot_prnp_across_samples(
    cellrouter = cellrouter_neftel,
    dataset = "neftel",
    plot_width = 9,
    plot_height = 3
)

## Step 3: Establishing groups with distinct PRNP levels -------------------

# Classifying cells as PRNP+ or PRNP-.
cellrouter_neftel <- add_prnp_status(cellrouter_neftel)

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.neftel <- concatenate_data(cellrouter = cellrouter_neftel)

#### @ FIGURE 2B, CENTRAL PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
plot_prnp_umap(
    dim.red = dim.red.neftel,
    dataset = "neftel",
    plot_width = 4.5,
    plot_height = 4
)

#### @ FIGURE S2A, CENTRAL PANEL (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample.
plot_sample_umap(
    dim.red = dim.red.neftel,
    dataset = "neftel",
    plot_width = 4.5,
    plot_height = 4
)

## Step 4: Finding signatures of PRNP+ and PRNP- cells ---------------------

markers_neftel <- get_signatures(
    cellrouter = cellrouter_neftel, 
    dataset = "neftel"
)

# Currently not used in the manuscript.
# plot_signature_volcano(markers_neftel, "neftel")

# Perform GSEA on markers found for PRNP positive cells.
gsea_neftel_prnp_pos <- do_gsea(
    df = markers_neftel,
    dataset = "neftel",
    population_ = "PRNP_positive_cells"
)
plot_gsea(
    gse = gsea_neftel_prnp_pos,
    dataset = "neftel",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform GSEA on markers found for PRNP negative cells.
gsea_neftel_prnp_neg <- do_gsea(
    df = markers_neftel,
    dataset = "neftel",
    population_ = "PRNP_negative_cells"
)
plot_gsea(
    gse = gsea_neftel_prnp_neg,
    dataset = "neftel",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP positive cells.
ora_neftel_prnp_pos <- do_ora(
    df = markers_neftel,
    dataset = "neftel",
    population_ = "PRNP_positive_cells"
)
plot_ora(
    gse = ora_neftel_prnp_pos,
    dataset = "neftel",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP negative cells.
ora_neftel_prnp_neg <- do_ora(
    df = markers_neftel,
    dataset = "neftel",
    population_ = "PRNP_negative_cells"
)
plot_ora(
    gse = ora_neftel_prnp_neg,
    dataset = "neftel",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

## Step 5: Correlation analysis --------------------------------------------

df_neftel <- get_correlations_with_prnp(
    cellrouter = cellrouter_neftel, 
    dataset = "neftel"
)

#### @ FIGURE 2C, CENTRAL PANEL (MAIN) @ #### 
plot_correlations_with_prnp(
    df = df_neftel, 
    dataset = "neftel",
    plot_width = 5,
    plot_height = 3
)


# Richards et al. dataset -------------------------------------------------


# Data downloaded from:
# https://singlecell.broadinstitute.org/single_cell/study/SCP503/gradient-of-developmental-and-injury-reponse-transcriptional-states-define-functional-vulnerabilities-underpinning-glioblastoma-heterogeneity


## Step 1: Loading and subsetting data for analysis ------------------------

GBM_44k_raw_data <- read.csv(
    paste0(getwd(), "/data/richards/Richards_NatureCancer_GBM_scRNAseq_counts.csv.gz"),
    header = TRUE,
    row.names = 1,
    sep = ",",
    check.names = FALSE
)

GSC_and_whole_tumor_metadata <- read.table(
    paste0(getwd(), "/data/richards/GSCs_Tumour_MetaData.txt"),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

rownames(GSC_and_whole_tumor_metadata) <- gsub('-', '.', rownames(GSC_and_whole_tumor_metadata))
rownames(GSC_and_whole_tumor_metadata) <- gsub('GBM_', '', rownames(GSC_and_whole_tumor_metadata))

# Subset metadata to keep only malignant cells.
malignant_metadata_richards <- GSC_and_whole_tumor_metadata[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR",]

# Create new sample column that matches the study's clinical metadata sample IDs.
malignant_metadata_richards$Sample_clean <- sapply(1:nrow(malignant_metadata_richards), function(x) {
    if(!grepl("-", malignant_metadata_richards$Sample.ID[x])) {
        return(malignant_metadata_richards$Sample.ID[x])
    } else {
        sample_clean <- gsub("[-]([A-Z]).*", "", malignant_metadata_richards$Sample.ID[x])
        sample_clean <- paste0(sample_clean, "_T")
        return(sample_clean)
    }    
})

# Filter out samples that should not be used for downstream analyses.
malignant_metadata_richards <- filter_out_samples(
    mdata = malignant_metadata_richards, 
    sample_column = "Sample_clean"
)

# Subset count data to keep only malignant cells & cells from samples that should be kept for
# downstream analyses.
malignant_data_richards <- GBM_44k_raw_data[,colnames(GBM_44k_raw_data) %in% rownames(GSC_and_whole_tumor_metadata)[GSC_and_whole_tumor_metadata$Sample.Type == "TUMOUR"]]
malignant_data_richards <- malignant_data_richards %>% dplyr::select(rownames(malignant_metadata_richards))

## Step 2: Data processing -------------------------------------------------

# Create CellRouter object for Richards et al dataset.
cellrouter_richards <- CreateCellRouter(
    as.matrix(malignant_data_richards), 
    min.genes = 0,
    min.cells = 0,
    is.expr = 0
)

# Process Richards' CellRouter object.
cellrouter_richards <- process_cellrouter(
    cellrouter = cellrouter_richards,
    metadata = malignant_metadata_richards,
    normalized = FALSE,
    sample_column_id = "Sample_clean",
    dataset = "richards"
)

#### @ FIGURE ? (SUPPLEMENTAL) @ #### 
# Boxplot showing the normalized PRNP expression in each sample of the dataset.
plot_prnp_across_samples(
    cellrouter = cellrouter_richards,
    dataset = "richards",
    plot_width = 2.5,
    plot_height = 3
)

## Step 3: Establishing groups with distinct PRNP levels -------------------

# Classifying cells as PRNP+ or PRNP-.
cellrouter_richards <- add_prnp_status(cellrouter_richards)

# Concatenate dataframes to have cell embeddings, gene expression, and metadata information in the same table.
dim.red.richards <- concatenate_data(cellrouter = cellrouter_richards)

#### @ FIGURE 2B, RIGHT PANEL (MAIN) @ #### 
# UMAP showing cell colored by PRNP expression.
plot_prnp_umap(
    dim.red = dim.red.richards,
    dataset = "richards",
    plot_width = 4.5,
    plot_height = 4
)

#### @ FIGURE S2A, RIGHT PANEL (SUPPLEMENTAL) @ #### 
# UMAP showing cell colored by sample.
plot_sample_umap(
    dim.red = dim.red.richards,
    dataset = "richards",
    plot_width = 4.5,
    plot_height = 4
)

## Step 4: Finding signatures of PRNP+ and PRNP- cells ---------------------

markers_richards <- get_signatures(cellrouter_richards, "richards")

# Currently not used in the manuscript.
# plot_signature_volcano(markers_richards, "richards")

# Perform GSEA on markers found for PRNP positive cells.
gsea_richards_prnp_pos <- do_gsea(
    df = markers_richards,
    dataset = "richards",
    population_ = "PRNP_positive_cells"
)
plot_gsea(
    gse = gsea_richards_prnp_pos,
    dataset = "richards",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform GSEA on markers found for PRNP negative cells.
gsea_richards_prnp_neg <- do_gsea(
    df = markers_richards,
    dataset = "richards",
    population_ = "PRNP_negative_cells"
)
plot_gsea(
    gse = gsea_richards_prnp_neg,
    dataset = "richards",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP positive cells.
ora_richards_prnp_pos <- do_ora(
    df = markers_richards,
    dataset = "richards",
    population_ = "PRNP_positive_cells"
)
plot_ora(
    gse = ora_richards_prnp_pos,
    dataset = "richards",
    population_ = "PRNP_positive_cells",
    title = expression("Markers of"~italic("PRNP")~"positive cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

# Perform ORA on markers found for PRNP negative cells.
ora_richards_prnp_neg <- do_ora(
    df = markers_richards,
    dataset = "richards",
    population_ = "PRNP_negative_cells"
)
plot_ora(
    gse = ora_richards_prnp_neg,
    dataset = "richards",
    population_ = "PRNP_negative_cells",
    title = expression("Markers of"~italic("PRNP")~"negative cells"),
    terms = ,
    plot_width = 10,
    plot_height = 8
)

## Step 5: Correlation analysis --------------------------------------------

df_richards <- get_correlations_with_prnp(cellrouter_richards, "richards")

#### @ FIGURE 2C, RIGHT PANEL (MAIN) @ #### 
plot_correlations_with_prnp(
    df = df_richards,
    dataset = "richards",
    plot_width = 5,
    plot_height = 3
)
 

# Comparing the datasets ---------------------------------------------------

## Step 1: Finding common upregulated marker genes -------------------------

# Create list with upregulated markers of PRNP positive cells of each dataset
markers_prnp_positive_upGenes <- list(
    "darmanis" = markers_darmanis$gene[markers_darmanis$population == "PRNP_positive_cells" & markers_darmanis$fc > 0],
    "neftel" = markers_neftel$gene[markers_neftel$population == "PRNP_positive_cells" & markers_neftel$fc > 0],
    "richards" = markers_richards$gene[markers_richards$population == "PRNP_positive_cells" & markers_richards$fc > 0]
)

# Create Venn diagram.
venn.diagram(
    x = markers_prnp_positive_upGenes,
    category.names = c("Darmanis" , "Neftel" , "Richards"),
    filename = 'common_markers_prnp_positive_upGenes.png',
    output = TRUE,
    imagetype = "png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = brewer.pal(3, "Pastel1")
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
)

# Find the intersection of all datasets.
common_prnp_positive_upGenes <- Reduce(intersect, markers_prnp_positive_upGenes)

## Step 2: Finding common positively correlated genes ----------------------

# Create list with genes positively correlated with PRNP of each dataset
genes_positiveCorr_with_prnp <- list(
    "darmanis" = df_darmanis$gene2[df_darmanis$status == "Positive correlation"],
    "neftel" = df_neftel$gene2[df_neftel$status == "Positive correlation"],
    "richards" = df_richards$gene2[df_richards$status == "Positive correlation"],
)

# Create Venn diagram.
venn.diagram(
    x = genes_positiveCorr_with_prnp,
    category.names = c("Darmanis" , "Neftel" , "Richards"),
    filename = 'common_genes_positiveCorr_with_PRNP.png',
    output = TRUE,
    imagetype = "png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = brewer.pal(3, "Pastel1")
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
)

# Find the intersection of all datasets.
common_genes_positiveCorr_with_prnp <- Reduce(intersect, genes_positiveCorr_with_prnp)

# The End
sessionInfo()