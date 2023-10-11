# plotDRExpression --------------------------------------------------------

# The fusca_1.3.0 is not the most recent FUSCA version, therefore some functions
# need to be overwritten to return the plots instead of just showing them.

#' Show gene expression in the space of reduced dimensionality.
#'
#' Plot selected gene expression in reduced dimension space and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to show.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#' @param threshold numeric; threshold to rescale gene expression.
#' @param dims.use numeric vector; dimensions to use.
#' @param columns numeric; number of columns in the output figure.
#' @param dotsize numeric; dot size.
#' @param alpha numeric; transparency (between 0 and 1).
#' @param title character; the title of the plot.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname plotDRExpression-methods
setGeneric("plotDRExpression", function(object, assay.type='RNA', genelist,
                                        reduction.type = c('tsne', 'pca', 'DC',
                                                           'umap', 'custom'),
                                        threshold = 2, dims.use = c(1, 2),
                                        columns = 5, dotsize = 1, alpha = 0.5,
                                        title)
  standardGeneric("plotDRExpression"))
#' @rdname plotDRExpression-methods
#' @aliases plotDRExpression
setMethod("plotDRExpression",
          signature="CellRouter",
          definition=function(object, assay.type, genelist,
                              reduction.type = c('tsne', 'pca', 'DC', 'umap',
                                                 'custom'),
                              threshold = 2, dims.use = c(1, 2), columns = 5,
                              dotsize, alpha, title){
            reduction.type <- match.arg(reduction.type)
            matrix <- as.data.frame(slot(object, reduction.type)$
                                      cell.embeddings[ , dims.use])
            plots <- list()
            scores <- matrix
            colnames(scores) <- c('Dim_1', 'Dim_2')
            if(reduction.type == 'tsne'){
              xlab <- paste('tSNE ', dims.use[1], sep=' ')
              ylab <- paste('tSNE ', dims.use[2], sep=' ')
            } else if (reduction.type == 'pca'){
              xlab <- paste('PC', dims.use[1], sep='')
              ylab <- paste('PC', dims.use[2], sep='')
            } else if (reduction.type == 'DC'){
              xlab <- paste('DC', dims.use[1], sep='')
              ylab <- paste('DC', dims.use[2], sep='')
            }else if(reduction.type == 'umap'){
              xlab <- paste('UMAP', dims.use[1], sep='')
              ylab <- paste('UMAP', dims.use[2], sep='')
            } else {
              xlab <- 'Dim 1'
              ylab <- 'Dim 2'
            }
            # x <- slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(matrix)]
            x <- center_with_threshold(
              slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(matrix),
                                                         drop=FALSE],
              threshold)
            dfs <- data.frame()
            for(gene in genelist){
              expr <- x[gene, ]
              scores$GENE <- as.numeric(expr)
              scores$gene <- gene
              dfs <- rbind(dfs, scores)
            }
            dfs <- dfs[order(dfs$GENE), ]
            dfs$gene <- factor(dfs$gene, levels = genelist)
            # Create plot.
            p1 <- ggplot2::ggplot(dfs, ggplot2::aes(x = Dim_1, y = Dim_2,
                                                    colour = GENE)) +
              ggplot2::geom_point(size = dotsize) + # ggplot2::aes(alpha = GENE), 
              ggplot2::theme_bw() +
              ggplot2::scale_colour_gradient2("Relative expression",
                                              low = "red",
                                              mid = "white",
                                              high = "blue",
                                              midpoint = 0) +
              ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
              ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = "white"),
                             strip.background = ggplot2::element_blank()) +
              ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank()) +
              ggplot2::theme(# legend.position="none",
                             strip.background = ggplot2::element_rect(colour = "white",
                                                                      fill = "white")) +
              ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top",
                                                                title.hjust = 0.5)) +
                              # size = ggplot2::guide_legend(title.position = "top",
                              #                              title.hjust = 0.5)) +
              ggplot2::facet_wrap(~gene, ncol = columns) +
              ggplot2::ggtitle(title)
            return(p1)
          }
)




# plotReducedDimension ----------------------------------------------------


#' Plot reduced dimension space.
#'
#' Plot reduced dimension space using pca, tsne, diffusion components, UMAP or
#' custom.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#' @param dims.use numeric vector; dimensions to use.
#' @param annotation character; column in the metadata table to annotate the figure.
#' @param annotation.color character; corresponding color column for the annotation.
#' @param showlabels logical; show labels in the plot.
#' @param dotsize numeric; dot size.
#' @param labelsize numeric; label size.
#' @param convex boolean; whether to connect the dots.
#' @param percentage numeric; proportion of data in the panel (between 0 and 1).
#' @param alpha numeric; transparency (between 0 and 1).
#'
#' @return ggplot2; plot.
#'
#' @export
plotReducedDimension <- function(object, assay.type='RNA',
                                 reduction.type = c('tsne', 'pca', 'DC',
                                                    'umap', 'custom'),
                                 dims.use = c(1, 2), annotation,
                                 annotation.color, showlabels,
                                 dotsize = 1,
                                 labelsize = 3, convex = FALSE,
                                 percentage = 0.80, alpha = 0.1){
  reduction.type <- match.arg(reduction.type)
  # Get the dimensions from the reduction type slot, cells as rows.
  matrix <- slot(object, reduction.type)$cell.embeddings[ , dims.use]
  scores <- as.data.frame(matrix)
  # Uses the rows from the metadata information.
  scores <- scores[rownames(slot(object, 'assays')[[assay.type]]@sampTab), ]
  colnames(scores) <- c('x', 'y')
  # Factors the information from the annotation at the metadata table.
  scores$group <- factor(as.vector(
    slot(object, 'assays')[[assay.type]]@sampTab[[annotation]]),
    levels = unique(as.vector(
      slot(object, 'assays')[[assay.type]]@sampTab[[annotation]])))
  # Select colors.
  ## colors <- unique(slot(object, 'assays')[[assay.type]]@sampTab[[annotation.color]])
  colors <- c('blue', 'red')
  names(colors) <- unique(as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[annotation]]))
  # Create plot
  g1 <- ggplot2::ggplot(scores, ggplot2::aes(x = x, y = y, colour = group)) +
    ggplot2::theme(legend.position = 'right') + ggplot2::geom_point(size = dotsize)
  if(showlabels == TRUE){
    # Gets the median of each group.
    centers <- scores %>% dplyr::group_by(group) %>%
      dplyr::summarize(x = median(x = x), y = median(x = y))
    # Adding text for the top 20 genes.
    g1 <- g1 +
      ggplot2::geom_point(data = centers, mapping = ggplot2::aes(x = x, y = y),
                          size = 0, alpha = 0) +
      ggrepel::geom_text_repel(data = centers,
                               mapping = ggplot2::aes(label = group),
                               size = labelsize, colour = 'black')
  }
  if(reduction.type == 'tsne'){
    xlab <- paste('tSNE ', dims.use[1], sep=' ')
    ylab <- paste('tSNE ', dims.use[2], sep=' ')
  }else if(reduction.type == 'pca'){
    xlab <- paste('PC', dims.use[1], sep='')
    ylab <- paste('PC', dims.use[2], sep='')
  }else if(reduction.type == 'dc'){
    xlab <- paste('DC', dims.use[1], sep='')
    ylab <- paste('DC', dims.use[2], sep='')
  }else if(reduction.type == 'umap'){
    xlab <- paste('UMAP', dims.use[1], sep='')
    ylab <- paste('UMAP', dims.use[2], sep='')
  }else{
    xlab <- 'Dim 1'
    ylab <- 'Dim 2'
  }
  # Add theme, labels...
  g1 <- g1 + ggplot2::theme_bw() +
    ggplot2::ggtitle('')  +  ggplot2::scale_color_manual("", values = colors) +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = "white"),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank()) +
    ggplot2::guides(
      col = ggplot2::guide_legend(direction = "vertical", keywidth = 0.75,
                                  keyheight = 0.85, override.aes = list(size = 3)))
  # If convex.
  if(convex){
    # Prop is the proportion of data in the pannel.
    g1 <- g1 + stat_bag(prop = percentage, alpha = alpha)
  }
  return(g1)
}




# center_with_threshold ---------------------------------------------------

#' Center the data using threshold.
#'
#' Helper function of plotSignaturesHeatmap, plotPathHeatmap, plotClusters...
#'
#' @param center_data data frame; data to be processed.
#' @param threshold numeric; max value present in center_data.
center_with_threshold <- function(center_data, threshold){
  # Center data (automatically ignores zeros).
  center_data <- center_data - Matrix::rowMeans(center_data, na.rm = TRUE)
  # Keep values between threshold and -threshold.
  center_data[center_data > threshold] <- threshold
  center_data[center_data < (-1 * threshold)] <- -1 * threshold
  return(center_data)
}