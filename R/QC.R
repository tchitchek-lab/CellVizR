#' @title Computes the percentage of cell clusters with low number of cells
#'
#' @description This function aims to compute and show cell clusters having a number of associated cells lower than a specific threshold
#'
#' @param Celldata a Celldata object
#' @param th.size a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster
#' @param plot.device a boolean value specifying a results representation must be displayed
#'
#' @return a numerical value corresponding to the percentage of cell cluster with low number of cells
#'
#' @export
#'
QCSmallClusters <- function(Celldata,
                            th.size = 50,
                            plot.device = TRUE) {

  checkmate::qassert(th.size, "N1")
  checkmate::qassert(plot.device, "B1")

  values.small <- computeSmallClusters(Celldata,
                                       th.size = th.size)

  if (plot.device == TRUE) {
    plot <- plotSmallClusters(values.small)
    plot(plot)
  }

  return(values.small)
}

# @title Internal - Computes the percentage of clusters with low number of cells
#
# @description This function is used internally to compute the percentage of clusters having a number of associated cells lower than a specific threshold
#
#
# @param Celldata a Celldata object
# @param th.size a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster
#
# @return a list containing QC information for small clusters
# Returns a data.frame with a boolean value indicating if the number of associated cells is greater or less than the threshold
# In addition, returns the percentage of the calculation
#
# @export
#
computeSmallClusters <- function(Celldata,
                                 th.size) {

  data.clusters <- Celldata@matrix.cell.count

  total.cells <- apply(data.clusters, 1,
                       FUN = function(x) {
                         ifelse((sum(x) < th.size), TRUE, FALSE) })
  data.clusters[data.clusters < th.size] <- TRUE
  data.clusters[data.clusters >= th.size] <- FALSE
  data.clusters <- cbind(data.clusters, total.cells = total.cells)
  data.clusters <- apply(data.clusters, 2, as.logical)
  perc <- sum(data.clusters[, "total.cells"] == TRUE) / nrow(data.clusters) * 100

  return(list(data.clusters = data.clusters, perc = perc))
}

# @title Internal - Plots a representation of QC for small clusters
#
# @description This function is used internally to create a representation showing the fraction of clusters having a number associated cell lower than a specific threshold
#
# @param values.small a list providing the small cluster QC information. Such as data.frame containing the boolean values and the percentage computed
#
# @return a ggplot2 object
#
plotSmallClusters <- function(values.small) {

  data.clusters <- values.small$data.clusters
  perc <- values.small$perc

  data.melted <- reshape2::melt(data.clusters)
  colnames(data.melted) <- c("clusters", "samples", "small")
  data.melted$clusters <- factor(data.melted$clusters,
                                 levels = rev(gtools::mixedsort(unique(data.melted$clusters))))

  title <- paste("small clusters quality control")
  subtitle <- paste0("percentage of clusters having a small number of cells = ",
                     format(round(perc, 2), nsmall = 2), "%")

  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
    ggplot2::geom_tile(data = data.melted,
                       ggplot2::aes_string(x = "samples", y = "clusters",
                                           fill = "small"), colour = "black") +
    ggplot2::scale_fill_manual(labels = c("TRUE" = "Small clusters"),
                               values = c("TRUE" = "red"),
                               na.value = "grey60") +
    ggplot2::geom_vline(xintercept = (length(unique(data.melted$samples)) - 0.5),
                        colour = "black", size = 2)

  plot <- plot +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::xlab("samples") +
    ggplot2::ylab("clusters")

  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.line = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank())

  return(plot)
}

#' @title Computes the percentage of clusters with uniform phenotypes
#'
#' @description This function aims to identify and show cell clusters having a uniform phenotype
#'
#' A uniform cluster corresponds to a cluster that have a unimodal expression and a low spread of expression for all its markers
#'
#' @details
#' -'uniform' corresponds to the verification of the unimodal distribution of markers with a Hartigans test
#'
#' -'IQR' corresponds to the verification of the distribution of markers so that they are not below the IQR threshold (interquantile range)
#'
#' -'both' corresponds to the combination of the two parameters : uniform and IQR
#'
#' @param Celldata a Celldata object
#' @param uniform.test a character providing the name of test assessment to perform. Possible value are : 'both', 'uniform', 'IQR'
#' @param th.pvalue a numeric value providing the p-value threshold of the Hartigan's dip test (unimodal if pvalue > th.pvalue)
#' @param th.IQR a numeric value providing the IQR (interquartile range) threshold to assume a distribution as uniform
#' @param plot.device a boolean value specifying if result representation must be displayed
#'
#' @return a numerical value corresponding to the percentage of cell cluster with unimodal expression and a low spread
#'
#' @export
#'
QCUniformClusters <- function(Celldata,
                              uniform.test = c("both", "uniform", "IQR"),
                              th.pvalue = 0.05,
                              th.IQR = 2,
                              plot.device = TRUE) {

  uniform.test <- match.arg(uniform.test)

  checkmate::qassert(uniform.test, "S1")
  checkmate::qassert(th.pvalue, "N1")
  checkmate::qassert(th.IQR, "N1")
  checkmate::qassert(plot.device, "B1")

  values.uniform <- computeUniformClusters(Celldata,
                                           uniform.test = uniform.test,
                                           th.pvalue = th.pvalue,
                                           th.IQR = th.IQR)

  if (plot.device == TRUE) {
    plot <- plotUniformClusters(values.uniform)
    plot(plot)
  }

  return(values.uniform)
}

# @title Internal - Computes percentage of clusters with uniform phenotype
#
# @description This function is used internally to identify cell clusters that have a non-uniform phenotype.
# A uniform cluster corresponds to clusters that have a unimodal expression and having a low spread of expression for all the markers to compose it
#
# @param Celldata a Celldata object
# @param uniform.test a character providing the name of test assessment to perform. Possible value are : 'both', 'uniform', 'IQR'
# @param th.pvalue a numeric value providing the p-value threshold of the Hartigan's dip test (unimodal if pvalue > th.pvalue)
# @param th.IQR a numeric value providing the IQR (interquartile range) threshold to assume a distribution as uniform
#
# @return a list containing QC information for small clusters
# Returns a data.frame with a boolean value indicating if the clusters have a uniform phenotype
# In addition, returns the percentage of the calculation
#
computeUniformClusters <- function(Celldata,
                                   uniform.test,
                                   th.pvalue,
                                   th.IQR) {

  computemode <- function(x) {
    den <- stats::density(x, kernel = c("gaussian"))
    return(list(x = den$x[den$y == max(den$y)], y = max(den$y)))
  }

  diptests <- data.frame()
  IQRs <- data.frame()
  modes <- data.frame()
  for (cluster in gtools::mixedsort(unique(Celldata@identify.clusters))) {
    for (marker in Celldata@identify.clusters.params$clustering.markers) {
      matrix.expression.sub <- Celldata@matrix.expression[Celldata@identify.clusters == cluster, marker]
      diptest <- diptest::dip.test(matrix.expression.sub)$p.value
      quantiles <- stats::quantile(matrix.expression.sub)
      IQR <- quantiles[4] - quantiles[2]
      mode <- computemode(matrix.expression.sub)$y
      diptests[cluster, marker] <- diptest
      IQRs[cluster, marker] <- IQR
      modes[cluster, marker] <- mode
    }
  }

  diptests$clusters <- rownames(diptests)
  diptests.m <- reshape::melt(diptests, id.vars = "clusters")
  colnames(diptests.m) <- c("clusters", "markers", "pv_dip")

  IQRs$clusters <- rownames(IQRs)
  IQRs.m <- reshape::melt(IQRs, id.vars = "clusters")
  colnames(IQRs.m) <- c("clusters", "markers", "IQR")

  res <- merge(diptests.m, IQRs.m)

  res$passed <- FALSE

  if (uniform.test == "uniform") {
    res$passed[res$pv_dip > th.pvalue] <- TRUE
  }
  if (uniform.test == "IQR") {
    res$passed[res$IQR < th.IQR] <- TRUE
  }
  if (uniform.test == "both") {
    res$passed[res$pv_dip > th.pvalue & res$IQR < th.IQR] <- TRUE
  }

  QC <- plyr::ddply(res, "clusters", function(x) {
    all(x$passed)
  })
  perc <- sum(QC$V1 == TRUE) / nrow(QC) * 100

  QC$markers <- "whole phenotype"
  QC$pv_dip <- NA
  QC$IQR <- NA
  QC$passed <- QC$V1
  QC$V1 <- NULL
  res <- rbind(res, QC)

  res$clusters <- factor(res$clusters,
                         levels = rev(gtools::mixedsort(unique(res$clusters))))

  return(list(res = res, perc = perc))
}

# @title Internal - Plots a representation of QC for uniform phenotype cluster
#
# @description This function is used internally to create a graphic representation
#
# @param values.uniform a list providing the uniform cluster QC information. Such as data.frame containing the boolean values and the percentage computed
#
# @return a ggplot2 object
#
plotUniformClusters <- function(values.uniform) {

  res <- values.uniform$res
  perc <- values.uniform$perc

  title <- paste("Uniform clusters quality control")
  subtitle <- paste0("percentage of clusters having a uniform phenotype = ",
                     format(round(perc, 2), nsmall = 2), "%")

  plot <- ggplot2::ggplot(data = res) +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
    ggplot2::geom_tile(ggplot2::aes_string(x = "markers", y = "clusters",
                                           fill = "passed"), colour = "black") +
    ggplot2::scale_fill_manual(labels = c("FALSE" = "No uniform"),
                               values = c("FALSE" = "red"),
                               na.value = "grey70")

  plot <- plot +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::xlab("clustering markers") +
    ggplot2::ylab("clusters")

  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.line = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank())

  return(plot)
}
