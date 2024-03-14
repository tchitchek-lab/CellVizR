#' @title Computes the percentage of cell clusters with low number of cells
#'
#' @description This function aims to compute and show cell clusters having a number of associated cells lower than a specific threshold.
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
# @description This function is used internally to compute the percentage of clusters having a number of associated cells lower than a specific threshold.
#
# @param Celldata a Celldata object
# @param th.size a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster
#
# @return a list containing QC information for small clusters.
# Returns a data.frame with a boolean value indicating if the number of associated cells is greater or less than the threshold.
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
# @description This function is used internally to create a representation showing the fraction of clusters having a number associated cell lower than a specific threshold.
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
                       ggplot2::aes(x = samples, y = clusters,
                                           fill = small), colour = "black") +
    ggplot2::scale_fill_manual(labels = c("TRUE" = "Small clusters"),
                               values = c("TRUE" = "red"),
                               na.value = "grey60") +
    ggplot2::geom_vline(xintercept = (length(unique(data.melted$samples)) - 0.5),
                        colour = "black", linewidth = 2)
  
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
#' @description This function aims to identify and show cell clusters having a uniform phenotype.
#'
#' A uniform cluster corresponds to a cluster that have a unimodal expression and a low spread of expression for all its markers.
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
# A uniform cluster corresponds to clusters that have a unimodal expression and having a low spread of expression for all the markers to compose it.
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
# @description This function is used internally to create a graphic representation.
#
# @param values.uniform a list providing the uniform cluster QC information. Such as data.frame containing the boolean values and the percentage computed.
#
# @return a ggplot2 object
#
plotUniformClusters <- function(values.uniform) {
  
  res <- values.uniform$res
  perc <- values.uniform$perc
  
  title <- paste("Uniform clusters quality control")
  subtitle <- paste0("percentage of clusters having a uniform phenotype = ",
                     format(round(perc, 2), nsmall = 2), "%")
  
  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
    ggplot2::geom_tile(data = res, 
                       ggplot2::aes(x = markers, y = clusters,
                                           fill = passed), colour = "black") +
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

#' @title QC for dimensionality reduction : Computes correlation between pairwise distances in the high-dimensional space and in the embedding
#'
#' @description CPD quantifies preservation of the global, or macroscropic structure.
#' A performant RD should give a boxplot with linear relationships betweens pairwise distance in the different data spaces
#'
#' @param Celldata a S4 object of class 'Celldata'
#' @param downsampling a numeric being the sample size used for computation
#' @param method a character being the method use to compute correlation. Possible values are : "pearson", "kendall", "spearman". Default value is "pearson"
#' @param plot.device a boolean value specifying if result representation must be displayed
#'
#' @return a list containing correlation coefficient and boxplot of pairwise distances
#'
#' @export
#' 
QCCorrelationManifold <- function(Celldata, 
                                  downsampling = 1000, 
                                  method = c("pearson", "kendall", "spearman"), 
                                  plot.device = TRUE) {
  
  method = match.arg(method)
  checkmate::qassert(method, "S1")
  checkmate::qassert(downsampling, "N1")
  
  values.correlation <- computeCorrelationManifold(Celldata, 
                                                   downsampling = downsampling, 
                                                   method = method)
  
  if (plot.device == TRUE) {
    plot <- plotCorrelationManifold(values.correlation)
    plot(plot)
  }
  
  return(values.correlation)
  
}

# @title QC for dimensionality reduction : Computes correlation between pairwise distances in the high-dimensional space and in the embedding
#
# @description CPD quantifies preservation of the global, or macroscropic structure.
# A performant RD should give a boxplot with linear relationships betweens pairwise distance in the different data spaces
#
# @param Celldata a S4 object of class 'Celldata'
# @param subsetsize a numeric being the sample size used for computation
# @param method a character being the method use to compute correlation. Possible values are : "pearson", "kendall", "spearman". Default value is "pearson"
#
# @return a list containing correlation coefficient 

computeCorrelationManifold <- function(Celldata, 
                                       downsampling, 
                                       method) {
  
  if (nrow(Celldata@manifold) == 0) { 
    stop("Error in 'compute.KNclass.error' function : KNclass QC require dimensionality step to be performed but
  'DimReduction' slot from 'Celldata' argument is empty") }
  
  sampled = sample(seq(1, nrow(Celldata@manifold)), downsampling)
  
  dist.DR = as.vector(stats::dist((Celldata@manifold[sampled,])))
  dist.origin = as.vector(stats::dist(Celldata@matrix.expression[sampled,]))
  
  corr = stats::cor(dist.DR, dist.origin, method=method)
  
  return(list(dist.DR = dist.DR, dist.origin = dist.origin, corr = corr))
}

# @title QC for dimensionality reduction : Computes correlation between pairwise distances in the high-dimensional space and in the embedding
#
# @description CPD quantifies preservation of the global, or macroscropic structure.
# A performant RD should give a boxplot with linear relationships betweens pairwise distance in the different data spaces
#
# @param Celldata a S4 object of class 'Celldata'
# @param subsetsize a numeric being the sample size used for computation
# @param method a character being the method use to compute correlation. Possible values are : "pearson", "kendall", "spearman". Default value is "pearson"
#
# @return a list containing correlation coefficient and boxplot of pairwise distances
#

plotCorrelationManifold <- function(values.correlation) {
  
  dist.DR <- values.correlation$dist.DR
  dist.origin <- values.correlation$dist.origin 
  
  matrix.dist = data.frame("DR" = dist.DR, "OR" = cut(dist.origin, 50))
  
  title <- paste("Pairwise distances in original space vs reduced space")
  subtitle <- paste0("correlation of pairwise distance = ",
                     format(round(values.correlation$corr, 2), nsmall = 2))
  
  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
    ggplot2::geom_boxplot(data = matrix.dist,
                          ggplot2::aes(x = OR, 
                                              y = DR),
                          fill="slateblue", alpha=0.2)
  plot <- plot +
    ggplot2::xlab("Distance cuts in original space") + 
    ggplot2::ylab("Distance in reduced space")
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust=0.5, size = 20, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20))
  
  return(plot)
}

#' @title QC for dimensionality reduction : Computes proportion of nearest neighbours preservation in reduced dimension in comparison to high dimension
#'
#' @description The fraction of k-nearest neighbours in the original high-dimensional data that are preserved as k-nearest neighbours in the embedding
#' It compute the average across all n points. KNN quantifies preservation of the local, or microscopic structure.
#'
#' @param Celldata a S4 object of class 'Celldata'
#' @param KNN a numeric being the number of nearest neighbors to compute
#' @param downsampling a numeric being the sample size used for computation
#'
#' @return a numeric being proportion
#'
#' @export

QCKnnManifold <- function(Celldata, 
                          KNN = 5, 
                          downsampling = 50 ) {
  
  checkmate::qassert(downsampling, "N1")
  checkmate::qassert(KNN, "N1")
  
  if (nrow(Celldata@manifold) == 0) { 
    stop("Error in 'compute.KNclass.error' function : KNclass QC require dimensionality step to be performed but
  'DimReduction' slot from 'Celldata' argument is empty") }
  
  sampled = sample(seq(1, nrow(Celldata@manifold)), downsampling) 
  
  data.or = Celldata@matrix.expression[sampled, ]
  data.dr = Celldata@manifold[sampled, ]
  
  neighborMatrix.or = RANN::nn2(data.or, data.or, KNN+1, 
                                searchtype = "standard")[[1]][,-1]
  neighborMatrix.dr = RANN::nn2(data.dr, data.dr, KNN+1, 
                                searchtype = "standard")[[1]][,-1]
  
  proportion.vector = sapply(1:nrow(data.or),
                             function(i){
                               same.nn = intersect(neighborMatrix.or[i,], 
                                                   neighborMatrix.dr[i,])
                               return(length(same.nn))
                             })
  res = mean(proportion.vector)/KNN
  
  return(res)
}

#' @title QC for dimensionality reduction : Computes proportion of nearest clusters preservation in reduced dimension in comparison to high dimension
#'
#' @description The fraction of k-nearest class (clusters) means in the original data that are preserved as k-nearest class means in the embedding.
#' it computes the class means only and averages across all classes. KNC quantifies preservation of the mesoscopic structure.
#'
#' @param Celldata a S4 object of class 'Celldata'
#' @param KNC a numeric being the number of nearest classes to compute
#' @param downsampling a numeric being the sample size used for computation
#'
#' @return a numeric being proportion
#'
#' @export
#'

QCKncManifold <- function(Celldata, 
                          KNC = 5, 
                          downsampling = 50){
  
  checkmate::qassert(downsampling, "N1")
  checkmate::qassert(KNC, "N1")
  
  if (nrow(Celldata@manifold) == 0) { 
    stop("Error in 'compute.KNclass.error' function : KNclass QC require dimensionality step to be performed but
  'DimReduction' slot from 'Celldata' argument is empty") }
  if (length(Celldata@identify.clusters) == 0) { 
    stop("Error in 'compute.KNclass.error' function : KNclass QC require clustering step to be performed but
  'Clustering' slot from 'Celldata' argument is empty") }
  
  sampled = sample(seq(1,nrow(Celldata@manifold)), downsampling)
  
  data.or = cbind.data.frame("clusters" = Celldata@identify.clusters[sampled], Celldata@matrix.expression[sampled, ])
  data.dr = cbind.data.frame("clusters" = Celldata@identify.clusters[sampled], Celldata@manifold[sampled, ])
  
  class.means.or = plyr::ddply(data.or,
                               "clusters",
                               function(x){ 
                                 return(colMeans(x[,-1])) })
  
  class.means.dr = plyr::ddply(data.dr,
                               "clusters",
                               function(x){ 
                                 return(colMeans(x[,-1])) })
  
  neighborMatrix.or = RANN::nn2(class.means.or, class.means.or, 
                                KNC+1, searchtype = "standard")[[1]][,-1]
  neighborMatrix.dr = RANN::nn2(class.means.dr, class.means.dr, 
                                KNC+1, searchtype = "standard")[[1]][,-1]
  
  proportion.vector = sapply(seq(1, nrow(class.means.or)),
                             function(i){
                               same.nc = intersect(neighborMatrix.or[i,], neighborMatrix.dr[i,])
                               return(length(same.nc))
                             })
  
  res = mean(proportion.vector)/KNC
  
  return(res)
}
