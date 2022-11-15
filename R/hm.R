# @title Internal - Rescales marker expression by quantile
#
# @description This function is used internally to rescale the marker expression by the quantile method
#
# @param exprs a data.frame containing the marker expressions for each cluster
# @param quant.low a numeric value providing the value of the first quantile
# @param quant.high a numeric value providing the value of the last quantile
#
# @return a data.frame containing quantile rescale marker expressions
#
rescaleMarkerExpressions <- function(exprs,
                                    quant.low = 0.05,
                                    quant.high = 0.95) {

  exprs <- plyr::ddply(exprs, "marker", function(df) {
    quantiles <- stats::quantile(df$value, probs = c(quant.low, quant.high))
    low <- quantiles[1]
    high <- quantiles[2]
    df$value[df$value < low] <- low
    df$value[df$value > high] <- high
    values <- scales::rescale(df$value, from = c(low, high), to = c(0, 4))

    return(data.frame(clusters = df$clusters, value = values))
  })

  return(exprs)
}

# @title Internal - Computes marker expression medians for each cluster
#
# @description This function is used internally to computes marker expression medians for each cluster.
#
# @param exprs a data.frame containing the marker expressions for each cell
#
# @return a data.frame containing the median of marker expressions for each cluster (clusters, markers, medians)
#
computeMarkerMedians <- function(exprs) {

    exprs <- plyr::ddply(exprs, c("clusters", "marker"), function(df) {
    med <- stats::median(df$value)
    return(med) })

  colnames(exprs) <- c("clusters", "marker", "med")
  exprs$clusters <- as.character(exprs$clusters)
  exprs$marker <- as.vector(exprs$marker)

  return(exprs)
}

#' @title Plots an heatmap of cell marker expressions
#'
#' @description This function aims to visualize the cell marker expressions for selected markers and clusters
#'
#' The mean of median marker expressions is computed for each cluster, and marker expressions displayed using a categorical heatmap (5 categories are defined by default)
#' The range expression of each cell marker is discretized into several categories between bounds of marker expressions
#' To hierarchical clustering, shown using dendrogramm, can be computed on both marker and cluster levels
#'
#' @param Celldata a Celldata object
#' @param markers a character vector providing the marker names to use. By default, all markers are used
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param method.hclust a character value providing the agglomeration method to be use. Possible values are: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid' (please refer to the function 'hclust' of the 'stats' package)
#' @param nb.cat a numeric specifying the number of categories to use
#' @param seed a numeric value providing the random seed to use during stochastic operations
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotHmExpressions <- function(Celldata,
                             markers = NULL,
                             clusters = NULL,
                             method.hclust = c("ward.D", "ward.D2", "single",
                                               "complete", "average", "mcquitty",
                                               "median", "centroid"),
                             nb.cat = 5,
                             seed = 42) {

  method.hclust <- match.arg(method.hclust)

  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(clusters, c("0", "S*", "I*"))
  checkmate::qassert(method.hclust, "S1")
  checkmate::qassert(nb.cat, "N1")
  checkmate::qassert(seed, "N1")

  rainbow.color.palette <- c("white", "yellow", "orange", "red", "red4")
  rainbow.color.palette <- grDevices::colorRampPalette(rainbow.color.palette)(nb.cat)

  if (!is.null(markers)) {
    exp.markers <- Celldata@matrix.expression[, markers]
  } else {
    exp.markers <- Celldata@matrix.expression[, colnames(Celldata@matrix.expression)]
  }

  exp.markers$clusters <- Celldata@identify.clusters

  if (!is.null(clusters)) {
    exp.markers <- exp.markers[exp.markers$clusters %in% clusters, ]
  }

  proj.melt <- reshape::melt(exp.markers, id = "clusters")
  colnames(proj.melt) <- c("clusters", "marker", "value")

  proj.melt <- rescaleMarkerExpressions(proj.melt)
  proj.melt <- computeMarkerMedians(proj.melt)

  row.data <- reshape::cast(proj.melt, marker ~ clusters, value = "med")
  row.data$marker <- NULL
  hclust.row <- stats::hclust(stats::dist(row.data), method = method.hclust)
  labels.row <- hclust.row$labels[hclust.row$order]
  proj.melt$marker <- factor(proj.melt$marker, levels = labels.row)

  col.data <- reshape::cast(proj.melt, clusters ~ marker, value = "med")
  col.data$clusters <- NULL
  hclust.col <- stats::hclust(stats::dist(col.data), method = method.hclust)

  do.call("set.seed", list(seed))
  dist.col <- stats::dist(col.data)
  isoMDS <- MASS::isoMDS(dist.col, k = 1, trace = FALSE)
  col.tsne <- isoMDS$points

  col.order <- rownames(col.tsne)[order(col.tsne)]
  hclust.col <- dendextend::rotate(hclust.col, col.order)
  labels.col <- hclust.col$labels[hclust.col$order]
  proj.melt$clusters <- factor(proj.melt$clusters, levels = labels.col)

  plot <- ggplot2::ggplot(proj.melt,
                          ggplot2::aes_string(x = "clusters", y = "marker")) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "med"), color = "black", size = 0.1) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))

  plot <- plot +
    ggplot2::scale_fill_gradientn(colours = rainbow.color.palette, na.value = "black") +
    ggplot2::labs(fill = "Median")

  plot <- plot +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank())

  dend.data <- ggdendro::dendro_data(stats::as.dendrogram(hclust.row), type = "rectangle")
  dendro.row <- ggplot2::ggplot(dend.data$segments) +
    ggplot2::geom_segment(ggplot2::aes_string(x = "x", y = "y",
                                              xend = "xend", yend = "yend")) +
    ggplot2::scale_x_continuous(expand = c(0.5 / length(unique(proj.melt$marker)),
                                           0.5 / length(unique(proj.melt$marker)))) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void()
  label.row <- ggplot2::ggplot() +
    ggplot2::geom_text(data = dend.data$labels, ggplot2::aes_string(x = "x", label = "label"),
                       size = 3, y = 0, hjust = 0, vjust = 0.5, angle = 0) +
    ggplot2::scale_x_continuous(expand = c(0.5 / length(unique(proj.melt$marker)),
                                           0.5 / length(unique(proj.melt$marker)))) +
    ggplot2::coord_flip() +
    ggplot2::theme_void()

  dend.data <- ggdendro::dendro_data(stats::as.dendrogram(hclust.col), type = "rectangle")
  dendro.col <- ggplot2::ggplot(dend.data$segments) +
    ggplot2::geom_segment(ggplot2::aes_string(x = "x", y = "y",
                                              xend = "xend", yend = "yend")) +
    ggplot2::scale_x_continuous(expand = c(0.5 / length(unique(proj.melt$clusters)),
                                           0.5 / length(unique(proj.melt$clusters)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))
  label.col <- ggplot2::ggplot() +
    ggplot2::geom_text(data = dend.data$labels, ggplot2::aes_string(x = "x", label = "label"),
                       size = 2, y = 1, hjust = 0, vjust = 0.5, angle = -90) +
    ggplot2::scale_x_continuous(expand = c(0.5 / length(unique(proj.melt$clusters)),
                                           0.5 / length(unique(proj.melt$clusters)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none")

  grob <- list()
  grob[["hm"]] <- plot + ggplot2::theme(legend.position = "none")
  grob[["dendro.col"]] <- dendro.col
  grob[["dendro.row"]] <- dendro.row
  grob[["label.col"]] <- label.col
  grob[["label.row"]] <- label.row
  grob[["legend"]] <- cowplot::get_legend(plot)

  hm.exp <- gridExtra::arrangeGrob(
    grobs = grob,
    layout_matrix = rbind(
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  )
  

  hm.exp$order.clusters <- labels.col
  hm.exp$order.markers <- labels.row

  return(hm.exp)
}

#' @title Plots an heatmap of cell cluster abundances
#'
#' @description This function aims to visualize the abundances of cell clusters using an heatmap representation
#'
#' In such heatmap each column corresponds a cell cluster and he row corresponds the different samples
#' The heatmap can be restricted to specific cell clusters and samples.
#' The levels of abundance of each sample in each cluster is represented using a color gradient scale
#' Abundance values can be centered and reduced.
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param saturation a numeric value providing the saturation threshold of cell cluster abundances
#' @param rescale a boolean specifying if cell cluster abundances must be centered and reduced
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotHmAbundances <- function(Celldata,
                            clusters = NULL,
                            samples = NULL,
                            saturation = 2.5,
                            rescale = FALSE) {

  checkmate::qassert(clusters, c("0", "S*", "I*"))
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(saturation, "N1")
  checkmate::qassert(rescale, "B1")

  abundances <- Celldata@matrix.abundance

  if (rescale == TRUE) {
    abundances <- t(scale(t(abundances)))
    abundances <- data.frame(abundances)
  } else {
    abundances <- data.frame(abundances)
  }

  if (!is.null(samples)) {
    abundances <- abundances[, names(abundances) %in% samples]
  }

  abundances$clusters <- rownames(abundances)

  if (!is.null(clusters)) {
    abundances <- abundances[abundances$clusters %in% clusters, ]
    abundances$clusters <- factor(abundances$clusters, levels = clusters)
  } else {
    abundances$clusters <- factor(abundances$clusters, levels = gtools::mixedsort(unique(abundances$clusters)))
  }

  abundances.melt <- reshape::melt(abundances, id = "clusters")
  colnames(abundances.melt) <- c("clusters", "samples", "value")
  abundances.melt$samples <- as.vector(abundances.melt$samples)

  abundances.melt <- plyr::ddply(abundances.melt, c("samples", "clusters"), function(x) {
    return(mean(x$value))
  })
  colnames(abundances.melt) <- c("samples", "clusters", "value")

  abundances.melt$samples <- factor(abundances.melt$samples, levels = rev(unique(abundances.melt$samples)))

  abundances.melt$value[abundances.melt$value > saturation] <- saturation
  abundances.melt$value[abundances.melt$value < (-saturation)] <- -saturation

  abundances.melt$label <- as.numeric(as.vector(abundances.melt$clusters))
  abundances.melt$label <- stringr::str_pad(abundances.melt$label, 3, pad = "0")

  plot <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = abundances.melt,
                       ggplot2::aes_string(x = "clusters", y = "samples",
                                           fill = "value"),
                       color = NA, size = 2)

    plot <- plot +
    ggplot2::labs(fill = "abundance") +
    ggplot2::scale_fill_gradientn(colours = c("green", "black", "red"),
                                  limits = c(-saturation, saturation),
                                  breaks = c(-saturation, saturation)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))

  plot <- plot +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank())

  label.row <- ggplot2::ggplot() +
    ggplot2::geom_text(data = abundances.melt,
                       ggplot2::aes_string(x = "samples", label = "samples"),
                       size = 3, y = 0, hjust = 0, vjust = 0.5, angle = 0) +
    ggplot2::scale_x_discrete(expand = c(0.5 / length(unique(abundances.melt$samples)),
                                         0.5 / length(unique(abundances.melt$samples)))) +
    ggplot2::coord_flip() +
    ggplot2::theme_void()

  label.col <- ggplot2::ggplot() +
    ggplot2::geom_text(data = abundances.melt,
                       ggplot2::aes_string(x = "clusters", label = "clusters"),
                       size = 2, y = 1, hjust = 0, vjust = 0.5, angle = -90) +
    ggplot2::scale_x_discrete(expand = c(0.5 / length(unique(abundances.melt$clusters)),
                                         0.5 / length(unique(abundances.melt$clusters)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none")

  grob <- list()
  grob[["hm.abundances"]] <- plot + ggplot2::theme(legend.position = "none")
  grob[["label.row"]] <- label.row
  grob[["label.col"]] <- label.col
  grob[["legend"]] <- cowplot::get_legend(plot)

  hm.abundance <- gridExtra::arrangeGrob(
    grobs = grob,
    layout_matrix = rbind(
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,NA,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA,NA),
      c(NA,NA,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  )

  return(hm.abundance)
}

#' @title Plots an heatmap of a statistical analysis results
#'
#' @description This function aims to visualize the results of differential cell clusters analysis
#'
#' This representation displays statistical information for each cell cluster for a given comparison of samples
#' Different statistics can be visualized, such as the p-value, the log2(fold-change), and effect size
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param statistics a character value providing the name of the statistic to display. Possible values are: 'pvalue' for p-value, 'lfc' for log2 fold change or 'eff' for effect size
#' @param saturation a numeric value providing the saturation value for statistics to display
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotHmStatistics <- function(Celldata,
                            clusters = NULL,
                            statistics = c("pvalue", "lfc", "effsize"),
                            saturation = 3) {

  statistics <- match.arg(statistics)

  checkmate::qassert(clusters, c("0", "S*", "I*"))
  checkmate::qassert(statistics, "S1")
  checkmate::qassert(saturation, "N1")

  stats <- Celldata@statistic
  stats$value <- stats[, statistics]

  stats$clusters <- as.character(stats$clusters)

  if (!is.null(clusters)) {
    stats <- stats[stats$clusters %in% clusters, ]
    stats$clusters <- factor(stats$clusters, levels = clusters)
  } else {
    stats$clusters <- factor(stats$clusters, levels = gtools::mixedsort(unique(stats$clusters)))
  }

  if (statistics == "pvalue") {
    stats$value <- -log(stats$value) / log(10)
    stats$value <- stats$value * sign(stats$lfc)
  }

  stats$value[stats$value > saturation] <- saturation
  stats$value[stats$value < (-saturation)] <- -saturation

  stats$comparison <- factor(stats$comparison, levels = rev(unique(stats$comparison)))

  stats$label <- as.numeric(as.vector(stats$clusters))
  stats$label <- stringr::str_pad(stats$label, 3, pad = "0")

  plot <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = stats,
                       ggplot2::aes_string(x = "clusters", y = "comparison",
                                           fill = "value"),
                       color = "black", size = 0.1)

  plot <- plot +
    ggplot2::labs(fill = statistics) +
    ggplot2::guides(color = "none") +
    ggplot2::scale_fill_gradientn(limits = c(-saturation, saturation), colours = c("orange", "gray", "blue"),
                                  breaks = c(-saturation, 0, saturation),
                                  na.value = "gray", labels = function(x) {
                                    format(round(x, 1), nsmall = 1) }) +
    ggplot2::scale_color_manual(values = c("FALSE" = "NA", "TRUE" = "black")) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))

  plot <- plot +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  label.row <- ggplot2::ggplot() +
    ggplot2::geom_text(data = stats,
                       ggplot2::aes_string(x = "comparison", label = "comparison"),
                       size = 3, y = 0, hjust = 0, vjust = 0.5, angle = 0) +
    ggplot2::scale_x_discrete(expand = c(0.5 / length(unique(stats$comparison)),
                                         0.5 / length(unique(stats$comparison)))) +
    ggplot2::coord_flip() +
    ggplot2::theme_void()

  label.col <- ggplot2::ggplot() +
    ggplot2::geom_text(data = stats,
                       ggplot2::aes_string(x = "clusters", label = "clusters"),
                       size = 2, y = 1, hjust = 0, vjust = 0.5, angle = -90) +
    ggplot2::scale_x_discrete(expand = c(0.5 / length(unique(stats$clusters)),
                                         0.5 / length(unique(stats$clusters)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none")

  grob <- list()
  grob[["hm.stats"]] <- plot + ggplot2::theme(legend.position = "none")
  grob[["label.row"]] <- label.row
  grob[["label.col"]] <- label.col
  grob[["legend"]] <- cowplot::get_legend(plot)

  hm.stats <- gridExtra::arrangeGrob(
    grobs = grob,
    layout_matrix = rbind(
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,NA),
      c(NA,NA,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA,NA),
      c(NA,NA,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA),
      c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  )

  return(hm.stats)
}

#' @title Plots a combined expression and statistic heatmaps
#'
#' @description This function aims to combine the expression and statistic heatmaps.
#'
#' @param HM1 a ggplot object containing the expression heatmap
#' @param HM2 a ggplot object containing the statistic heatmap
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCombineHM <- function(HM1,
                          HM2) {

  grob2 <- c(HM1$grobs, HM2$grobs)

  combineHM <- gridExtra::grid.arrange(
    grobs = grob2,
    layout_matrix = rbind(
      c(6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA),
      c(6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA),
      c(6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5),
      c(NA,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,NA),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(11,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8),
      c(NA,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,NA),
      c(NA,NA,10,10,10,10,10,10,10,10,10,10,10,10,10,10,NA,NA))
  )

  return(combineHM)
}
