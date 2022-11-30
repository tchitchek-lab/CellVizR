#' @title Plots the numbers of cells for each sample
#'
#' @description This function aims to visualize the number of cells associated to each sample.
#'
#' This representation displays the samples in the X-axis and the number of associated cells in the Y-axis.
#' Several statistics can be computed and shown.
#'
#' @details
#' The following statistic can be computed:
#'
#' -'min' corresponds to the lowest number of cells within a data set
#'
#' -'median' corresponds to the number of cells separates the lower half from the upper half within data set
#'
#' -'mean' corresponds to the number of cells quantity shared within data set
#'
#' -'q75' corresponds to the number of cells separates the quantiles 75% within data set
#'
#' -'max' corresponds to the largest number of cells within a data set
#'
#' @param Celldata a Celldata object
#' @param stats a character vector providing the statistics to display. Possible values are: 'min', 'median', 'mean', 'q75', 'max'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param sort a boolean value indicating if clusters must be sorted by the number associated sample
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCellCounts <- function(Celldata,
                           stats = c("min", "median", "mean", "q75", "max"),
                           samples = NULL,
                           sort = TRUE) {
  
  checkmate::qassert(stats, "S*")
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(sort, "B1")
  
  nb.cells <- Celldata@matrix.expression
  nb.cells$samples <- Celldata@samples
  
  if (!is.null(samples)) {
    nb.cells <- nb.cells[nb.cells$samples %in% samples, ]
  }
  
  nb.cells <- plyr::ddply(nb.cells, "samples", nrow)
  colnames(nb.cells) <- c("samples", "nb")
  
  if (sort == TRUE) {
    nb.cells <- nb.cells[order(nb.cells$nb, decreasing = TRUE), ]
    nb.cells$samples <- factor(nb.cells$samples, levels = nb.cells$samples)
  } else {
    nb.cells$samples <- factor(nb.cells$samples, levels = nb.cells$samples)
  }
  
  min <- min(nb.cells$nb)
  mean <- mean(nb.cells$nb)
  median <- stats::median(nb.cells$nb)
  q75 <- stats::quantile(nb.cells$nb, probs = 0.75)
  max <- max(nb.cells$nb)
  
  xscale.middle <- nrow(nb.cells) / 2
  if (xscale.middle %% 2 == 0) {
    xscale.middle <- xscale.middle + 0.5
  }
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = nb.cells,
                          ggplot2::aes_string(x = "samples",
                                              xend = "samples",
                                              yend = "nb"),
                          y = 0, color = "gray50") +
    ggiraph::geom_point_interactive(data = nb.cells,
                                    ggplot2::aes_string(x = "samples",
                                                        y = "nb",
                                                        tooltip = "nb"),
                                    size = 2)
  
  if ("min" %in% stats) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = min, color = "blue") +
      ggplot2::geom_text(data = nb.cells,
                         ggplot2::aes_string(x = "xscale.middle",
                                             y = "min"),
                         label = paste0("min: ", format(min, scientific = FALSE,
                                                        big.mark = ",")),
                         color = "blue", vjust = 0)
  }
  
  if ("max" %in% stats) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = max, color = "blue") +
      ggplot2::geom_text(data = nb.cells,
                         ggplot2::aes_string(x = "xscale.middle",
                                             y = "max"),
                         label = paste0("max: ", format(max, scientific = FALSE,
                                                        big.mark = ",")),
                         color = "blue", vjust = 0)
  }
  
  if ("median" %in% stats) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = median, color = "red") +
      ggplot2::geom_text(data = nb.cells,
                         ggplot2::aes_string(x = "xscale.middle",
                                             y = "median"),
                         label = paste0("median: ", format(median, scientific = FALSE,
                                                           big.mark = ",")),
                         color = "red", vjust = 0)
  }
  
  if ("mean" %in% stats) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = mean, color = "orange") +
      ggplot2::geom_text(data = nb.cells,
                         ggplot2::aes_string(x = "xscale.middle",
                                             y = "mean"),
                         label = paste0("mean: ", format(mean, scientific = FALSE,
                                                         big.mark = ",")),
                         color = "orange", vjust = 0)
  }
  
  if ("q75" %in% stats) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = q75, color = "purple") +
      ggplot2::geom_text(data = nb.cells,
                         ggplot2::aes_string(x = "xscale.middle",
                                             y = "q75"),
                         label = paste0("q75: ", format(q75, scientific = FALSE,
                                                        big.mark = ",")),
                         color = "purple", vjust = 0)
  }
  
  plot <- plot +
    ggplot2::scale_y_continuous(labels = function(x) {
      format(x, scientific = FALSE,
             big.mark = ",")}) +
    ggplot2::ylab("number of cells")
  
  plot <- plot +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
                   axis.title.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "none")
  
  return(plot)
}

#' @title Plots the numbers of cells of each cluster
#'
#' @description This function aims to visualize the number of cells associated to each cluster.
#'
#' This representation displays the clusters in the X-axis and the total number of associated cells in the Y-axis.
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param sort a boolean value indicating if clusters must be sorted by the number associated cluster
#' @param legend.max.samples a numerical value specifying the maximal number of samples to display in the graphical legend
#'
#' @return a ggplot2 object
#' 
#' @export
#'
plotClustersCounts <- function(Celldata,
                               clusters = NULL,
                               sort = TRUE,
                               legend.max.samples = 10) {
  
  checkmate::qassert(clusters, c("0", "N+"))
  checkmate::qassert(sort, "B1")
  checkmate::qassert(legend.max.samples, "N1")
  
  matrix.cell <- Celldata@matrix.cell.count
  
  if (!is.null(clusters)) {
    matrix.cell <- matrix.cell[rownames(matrix.cell) %in% clusters, ]
  } else {
    clusters <- unique(Celldata@identify.clusters)
    matrix.cell <- matrix.cell[clusters, , drop = FALSE]
  }
  
  cells.number <- sum(colSums(matrix.cell))
  matrix.cell <- cbind(matrix.cell, "sum.of.samples" = apply(matrix.cell, 1, sum))
  
  matrix.cell <- cbind("clusters" = rownames(matrix.cell), matrix.cell)
  
  if (sort == TRUE) {
    order <- order(matrix.cell$sum.of.samples, decreasing = TRUE)
    matrix.cell$clusters <- factor(matrix.cell$clusters, levels = matrix.cell$clusters[order])
  } else {
    matrix.cell$clusters <- factor(matrix.cell$clusters, levels = clusters)
  }
  
  data.melted <- reshape2::melt(matrix.cell, id = c("clusters"))
  colnames(data.melted) <- c("clusters", "samples", "values")
  
  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(paste0("Count viewer (", format(cells.number,
                                                     big.mark = " "), " cells)", sep = "")) +
    ggiraph::geom_bar_interactive(data = data.melted[data.melted$samples != "sum.of.samples", ],
                                  ggplot2::aes_string(x = "clusters",
                                                      y = "values",
                                                      fill = "samples",
                                                      tooltip = "samples",
                                                      data_id = "samples"),
                                  stat = "identity", position = "stack", color = "black")
  
  plot <- plot +
    viridis::scale_fill_viridis(option = "turbo",
                                discrete = TRUE) +
    ggplot2::xlab("clusters") +
    ggplot2::ylab("number of cells")
  
  if (length(unique(Celldata@samples)) >= legend.max.samples) {
    plot <- plot + ggplot2::guides(fill = "none")
  }
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.key = ggplot2::element_blank())
  
  return(plot)
}

#' @title Plots of phenotype of identified cell clusters
#'
#' @description This function aims to visualize the density of marker/gene expressions for a set of given clusters.
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifier of the cluster to use
#' @param quant.low a numeric value providing the value of first quantile in the color gradiant
#' @param quant.high a numeric value providing the value of last quantile in the color gradiant
#' @param dip.th a numeric value specifying the p-value threshold of the Hartigan's dip test to consider non-unimodal clusters
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotMarkerDensity <- function(Celldata,
                              clusters,
                              quant.low = 0.05,
                              quant.high = 0.95,
                              dip.th = 0.01) {
  
  checkmate::qassert(clusters, "S+")
  checkmate::qassert(quant.low, "N1")
  checkmate::qassert(quant.high, "N1")
  checkmate::qassert(dip.th, "N1")
  
  expression.color.palette <- c("white", "yellow", "orange", "red", "red4")
  
  matrix.exp <- Celldata@matrix.expression
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(matrix.exp, cluster)
  
  parameters <- colnames(proj)
  markers <- parameters[!parameters %in% c("cluster")]
  
  grob <- list()
  
  for (marker in markers) {
    exp.values <- proj[, c(marker, "cluster")]
    colnames(exp.values) <- c("exp", "clusters")
    quantiles <- stats::quantile(exp.values$exp, probs = c(quant.low, quant.high))
    exp.values.clusters <- exp.values[exp.values$clusters %in% clusters, ]
    median <- median(exp.values.clusters$exp)
    tmp.exp <- exp.values.clusters$exp
    dip <- diptest::dip.test(tmp.exp)
    
    if (dip$p.value < dip.th) {
      color <- "red"
    } else {
      color <- "green"
    }
    
    plot <- ggplot2::ggplot() +
      ggplot2::ggtitle(marker) +
      ggridges::geom_density_ridges_gradient(data = exp.values,
                                             ggplot2::aes_string(x = "exp",
                                                                 y = "0",
                                                                 fill = "stat(x)")) +
      ggplot2::geom_density(data = exp.values.clusters,
                            ggplot2::aes_string(x = "exp"),
                            color = color,
                            fill = NA,
                            size = 0.75) +
      ggplot2::geom_vline(xintercept = median, linetype = "dashed", color = "gray20")
    
    plot <- plot +
      ggplot2::scale_fill_gradientn(limits = c(quantiles[1], quantiles[2]),
                                    colours = expression.color.palette,
                                    na.value = "black") +
      ggplot2::scale_x_continuous(limits = c(-2, 6), breaks = c(-2:6)) +
      ggplot2::scale_y_continuous(labels = function(x) {format(round(x, 1),nsmall = 1)})
    
    plot <- plot +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        legend.position = "none",
        legend.key.width = grid::unit(0.35, "cm"))
    
    grob[[marker]] <- plot
  }
  
  all.plot <- gridExtra::grid.arrange(grobs = grob)
  
  return(grid::grid.draw(all.plot))
}

#' @title Plots a representation of a computed manifold
#'
#' @description This function aims to visualize a computed manifold representation for given analysis.
#'
#' This representation can be used on a Celldata object for which a manifold analysis has been performed.
#'
#' If a cell clustering has been performed, then the clusters are delineated using concave hulls.
#' Additionally, the manifold can be colored based on the local cell density or marker expressions.
#' It is possible to centre and reduce the values of expressions.
#'
#' @param Celldata a Celldata object
#' @param markers a character value providing the name of the marker to use for the colouring. By default, cells are colored based on their local density
#' @param samples a character vector containing the names of biological samples to use
#' @param scale a boolean value specifying  if expression values must be rescaled
#' @param quant.low a numeric value providing the number of first quantile
#' @param quant.high a numeric value providing the number of last quantile 
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotManifold <- function(Celldata,
                         markers = "density",
                         samples = NULL,
                         scale = FALSE,
                         quant.low = 0.05,
                         quant.high = 0.95) {
  
  checkmate::qassert(markers, "S1")
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(quant.low, "N1")
  checkmate::qassert(quant.high, "N1")
  
  proj <- Celldata@manifold
  
  if (length(proj) == 0) {
    stop("Manifold is null. Please create manifold.")
  }
  
  colnames(proj) <- c("dim1", "dim2")
  clusters <- Celldata@identify.clusters
  
  if (length(clusters) != 0) {
    proj$clusters <- clusters
  }
  
  proj <- stats::na.omit(proj)
  lines <- Celldata@concave.hulls
  
  if (markers == "density") {
    legend.title <- "cell density"
    proj$value <- computeCellDensities(proj, n = 100)
  } else if (markers == "clusters") {
    legend.title <- "cell clusters"
    proj$value <- proj$clusters
  } else if (markers == "clusters" && length(clusters) == 0) {
    stop("Clusters is null")
  } else {
    legend.title <- "marker\nexpression"
    proj$value <- Celldata@matrix.expression[, markers]
  }
  
  plot <- ggplot2::ggplot()
  
  if (!is.null(samples)) {
    proj.ref <- proj
    proj <- proj[Celldata@samples %in% samples, ]
    plot <- plot +
      ggplot2::geom_point(data = proj.ref,
                          ggplot2::aes_string(x = "dim1",
                                              y = "dim2"),
                          color = "gray", size = 0.0001) +
      ggnewscale::new_scale_color()
  }
  
  plot <- plot +
    ggplot2::geom_point(data = proj,
                        ggplot2::aes_string(x = "dim1",
                                            y = "dim2",
                                            color = "value"),
                        size = 0.0001)
  
  if (nrow(lines) != 0) {
    plot <- plot + ggplot2::geom_path(data = lines,
                                      ggplot2::aes_string(x = "dim1",
                                                          y = "dim2",
                                                          group = "clusters"),
                                      color = "black", size = 0.2)
  }
  
  if (markers == "density") {
    plot <- plot + ggplot2::scale_color_gradientn(colours = c("black", "red"),
                                                  labels = function(x) {
                                                    return(rep("", length(x)))}) +
      ggplot2::labs(title = "density representation", color = legend.title)
    
  } else if (markers == "clusters") {
    
    proj.center <- plyr::ddply(proj, "clusters",function(x) {
      c(dim1 = stats::median(x$dim1),
        dim2 = stats::median(x$dim2))})
    
    plot <- plot + ggplot2::geom_text(data = proj.center,
                                      ggplot2::aes_string(x = "dim1",
                                                          y = "dim2",
                                                          label = "clusters"),
                                      size = 3) +
      # ggplot2::guides(color = "none") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=3))) +
      ggplot2::labs(title = "clusters representation", color = legend.title)
    
  } else {
    
    limits <- abs.proj(proj, scale, quant.low, quant.high)
    plot <- plot + viridis::scale_color_viridis(option = "plasma",
                                                direction = -1,
                                                limits = limits,
                                                oob = scales::oob_squish) +
      ggplot2::labs(title = markers, color = legend.title)
  }
  
  plot <- plot +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01))
  
  plot <- plot +
    ggplot2::xlab(paste0(Celldata@manifold.params$type, "1")) + 
    ggplot2::ylab(paste0(Celldata@manifold.params$type, "2")) + 
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   aspect.ratio = 1,
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(hjust = 1, vjust = 1, face = "bold"))
  
  return(plot)
}

#' @title Plots a PCA representation based cell cluster abundances
#'
#' @description This function aims to represent a Principal Component Analysis representation based on cell cluster abundances.
#' In such representation, clusters or samples are positioned based on computed principal components.
#' The representation can be displayed based on specific principal components.
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples.
#'
#' @param Celldata a Celldata object
#' @param levels a character value containing the variable to be displayed. Possible values are: 'both', 'clusters' or 'samples'
#' @param clusters a character vector containing the identifier of the cluster to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param components a numeric vector providing the components to display
#' @param condition.samples a character vector containing the variable to be studied for the samples. Possible values are: 'condition' or 'timepoint"
#' @param cor.radius.th a numeric value specifying the radius of the correlation plot radius
#' @param plot.text a boolean value specifying if adds text directly at the plot
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotPCA <- function(Celldata,
                    levels = c("both", "clusters", "samples"),
                    clusters = NULL,
                    samples = NULL,
                    components = c(1, 2),
                    condition.samples = c("condition", "timepoint"),
                    cor.radius.th = 0.6,
                    plot.text = TRUE) {
  
  levels <- match.arg(levels)
  condition.samples <- match.arg(condition.samples)
  
  checkmate::qassert(levels, "S1")
  checkmate::qassert(clusters, c("0", "S*"))
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(components, "N2")
  checkmate::qassert(condition.samples, "S1")
  checkmate::qassert(cor.radius.th, "N1")
  checkmate::qassert(plot.text, "B1")
  
  if (levels != "both" && levels != "clusters" && levels != "samples") {
    stop("The levels name is invalid")
  }
  
  circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
    r <- diameter / 2
    t <- seq(0, 2 * pi, length.out = npoints)
    x <- center[1] + r * cos(t)
    y <- center[2] + r * sin(t)
    return(data.frame(x = x, y = y))
  }
  
  matrix.abundance <- Celldata@matrix.abundance
  
  if (!is.null(clusters)) {
    matrix.abundance <- matrix.abundance[rownames(matrix.abundance) %in% clusters, ]
  }
  if (!is.null(samples)) {
    matrix.abundance <- matrix.abundance[, colnames(matrix.abundance) %in% samples]
  }
  
  res.PCA <- FactoMineR::PCA(t(matrix.abundance), graph = FALSE)
  var.explained <- res.PCA$eig[, 2]
  
  data.clusters <- data.frame(res.PCA$var$coord[, components])
  colnames(data.clusters) <- c("x", "y")
  data.clusters$id <- rownames(matrix.abundance)
  
  data.variables <- data.frame(res.PCA$ind$coord[, components])
  colnames(data.variables) <- c("x", "y")
  data.variables$id <- colnames(matrix.abundance)
  
  data.variables$x <- scales::rescale(data.variables$x, to = c(-1, 1))
  data.variables$y <- scales::rescale(data.variables$y, to = c(-1, 1))
  
  data.variables <- merge(data.variables, Celldata@metadata, by = "row.names")
  data.variables$Row.names <- NULL
  
  circle1 <- circleFun(c(0, 0), 2, npoints = 100)
  circle2 <- circleFun(c(0, 0), 2 * cor.radius.th, npoints = 100)
  
  data.clusters <- data.clusters[sqrt(data.clusters$x^2 + data.clusters$y^2) > cor.radius.th, ]
  
  xlab <- paste0("PCA", components[1], " (", round(var.explained[components[1]], 2), "%)")
  ylab <- paste0("PCA", components[2], " (", round(var.explained[components[2]], 2), "%)")
  
  plot <- ggplot2::ggplot()
  
  if (levels == "clusters" || levels == "both") {
    plot <- plot +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_path(data = circle1,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, color = "gray") +
      ggplot2::geom_path(data = circle2,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, linetype = "dashed", color = "gray") +
      ggiraph::geom_point_interactive(data = data.clusters,
                                      ggplot2::aes_string(x = "x", y = "y",
                                                          tooltip = "id",
                                                          data_id = "id"),
                                      shape = 21, size = 2,
                                      stroke = 0.1, fill = "black", color = "black")
    
    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = data.clusters,
                                 ggplot2::aes_string(x = "x",
                                                     y = "y",
                                                     label = "id"),
                                 size = 3,
                                 color = "black",
                                 max.overlaps = Inf,
                                 min.segment.length = 0,
                                 segment.color = NA,
                                 segment.size = 0.1)
    }
  }
  
  if (levels == "samples" || levels == "both") {
    plot <- plot +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_path(data = circle1,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, color = "gray") +
      ggplot2::geom_path(data = circle2,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, linetype = "dashed", color = "gray") +
      ggiraph::geom_point_interactive(data = data.variables,
                                      ggplot2::aes_string(x = "x", y = "y",
                                                          fill = condition.samples,
                                                          tooltip = "id",
                                                          data_id = "individual"),
                                      shape = 21, size = 2,
                                      stroke = 0.1, color = "black")
    
    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = data.variables,
                                 ggplot2::aes_string(x = "x", y = "y",
                                                     label = "id"),
                                 size = 3,
                                 color = "black",
                                 max.overlaps = Inf,
                                 min.segment.length = 0,
                                 segment.color = NA,
                                 segment.size = 0.1)
    }
  }
  
  plot <- plot + ggplot2::labs(title = "PCA representation")
  
  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank())
  
  return(plot)
}

#' @title Plots a MDS representation based on cell cluster abundances
#'
#' @description This function aims to visualize the similarities between samples or clusters based on their abundances, using a Multidimensional Scaling representation.
#' Each dot represents a sample or a cluster and the distances between the dots are proportional to the Euclidean distance between these objects.
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples.
#'
#' @param Celldata a Celldata object
#' @param matrix a character vector containing the matrix to be studied. Possible values are: 'abundance' or 'expression'
#' @param levels a character value containing the variable to be displayed. Possible values are: 'clusters' or 'samples'
#' @param condition.samples a character vector containing the variable to be studied for the samples. Possible values are: 'condition' or 'timepoint"
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param plot.text a boolean value specifying if adds text directly at the plot
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotMDS <- function(Celldata,
                    matrix = c("abundance","expression"),
                    levels = c("clusters", "samples"),
                    condition.samples = c("condition", "timepoint"),
                    clusters = NULL,
                    samples = NULL,
                    plot.text = TRUE) {
  
  levels <- match.arg(levels)
  condition.samples <- match.arg(condition.samples)
  matrix <- match.arg(matrix)
  
  checkmate::qassert(levels, "S1")
  checkmate::qassert(condition.samples, "S1")
  checkmate::qassert(matrix, "S1")
  checkmate::qassert(clusters, c("0", "S+"))
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(plot.text, "B1")
  
  if (levels != "clusters" && levels != "samples") {
    stop("The levels name is invalid")
  }
  
  if (matrix != "abundance" && matrix != "expression") {
    stop("The matrix name is invalid")
  }
  
  if (matrix == "abundance") {
    data.matrix <- Celldata@matrix.abundance
    if (!is.null(clusters)) {
      data.matrix <- data.matrix[rownames(data.matrix) %in% clusters, ]
    }
    if (!is.null(samples)) {
      data.matrix <- data.matrix[, names(data.matrix) %in% samples]
    }
  } else {
    data.matrix <- cbind(Celldata@matrix.expression, 
                         "samples" = Celldata@samples,
                         "clusters" = Celldata@identify.clusters)
    if (!is.null(clusters)){
      data.matrix = data.matrix[data.matrix$clusters %in% clusters,]
    }
    if (!is.null(samples)) {
      data.matrix = data.matrix[data.matrix$samples %in% samples,]
    }
    if (levels == "clusters") {
      data.matrix$samples <- NULL
      data.matrix <- plyr::ddply(data.matrix, "clusters", function(x) {
        x$clusters <- NULL
        apply(x, 2, stats::median)
      })
      rownames(data.matrix) <- data.matrix$clusters
      data.matrix$clusters <- NULL
    } else {
      data.matrix$clusters <- NULL
      data.matrix <- plyr::ddply(data.matrix, "samples", function(x) {
        x$samples <- NULL
        apply(x, 2, stats::median)
      })
      
      rownames(data.matrix) <- data.matrix$samples
      data.matrix$samples <- NULL
      data.matrix <- t(data.matrix)
    }
  }
  
  if (levels == "clusters") {
    
    dist.clusters <- stats::dist(data.matrix)
    fit1 <- MASS::isoMDS(dist.clusters, k = 2, trace = FALSE)
    
    x1 <- fit1$points[, 1]
    y1 <- fit1$points[, 2]
    
    min.lim <- min(min(x1), min(y1)) * 1.1
    max.lim <- max(max(x1), max(y1)) * 1.1
    
    proj.clusters <- data.frame(x = x1,
                                y = y1)
    
    proj.clusters$clusters <- rownames(data.matrix)
    
    plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ", round(fit1$stress, 2))) +
      ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggiraph::geom_point_interactive(data = proj.clusters,
                                      ggplot2::aes_string(x = "x", y = "y",
                                                          fill = "clusters",
                                                          tooltip = "clusters",
                                                          data_id = "clusters"),
                                      shape = 21, fill = "black", color = "black")
    
    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = proj.clusters,
                                 ggplot2::aes_string(x = "x",
                                                     y = "y",
                                                     label = "clusters"),
                                 size = 3)
    }
    
  } else {
    
    dist.samples <- stats::dist(t(data.matrix))
    
    fit2 <- MASS::isoMDS(dist.samples, k = 2, trace = FALSE)
    x2 <- fit2$points[, 1]
    y2 <- fit2$points[, 2]
    
    min.lim <- min(min(x2), min(y2)) * 1.1
    max.lim <- max(max(x2), max(y2)) * 1.1
    
    proj.samples <- data.frame(x = x2,
                               y = y2)
    
    proj.samples$samples <- colnames(data.matrix)
    
    proj.samples <- merge(proj.samples, Celldata@metadata, by = "row.names")
    proj.samples$Row.names <- NULL
    
    plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ", round(fit2$stress, 2))) +
      ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggiraph::geom_point_interactive(data = proj.samples,
                                      ggplot2::aes_string(x = "x",
                                                          y = "y",
                                                          fill = condition.samples,
                                                          tooltip = "samples",
                                                          data_id = "individual"),
                                      shape = 21, color = "black")
    
    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = proj.samples,
                                 ggplot2::aes_string(x = "x",
                                                     y = "y",
                                                     label = "samples"),
                                 size = 3)
    }
  }
  
  plot <- plot + ggplot2::coord_fixed()
  
  plot <- plot +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  return(plot)
}

#' @title Plots cell cluster abundances using a boxplot representation
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#'
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed.
#' The representation can be restricted to a specific set of samples.
#' Moreover, boxplot can be constructed based on sample meta information.
#' Statistic can be computed for all comparisons.
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to use
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observation a character value containing the parameters to use
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#' @param hide.ns a boolean value indicating if non-significant p-value must be hidden
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotBoxplot <- function(Celldata,
                        clusters,
                        samples = NULL,
                        observation = c("individual", "condition", "timepoint"),
                        value.y = c("percentage", "absolute"),
                        test.statistics = c("wilcox.test", "t.test"),
                        paired = FALSE,
                        hide.ns = TRUE) {
  
  observation <- match.arg(observation)
  test.statistics <- match.arg(test.statistics)
  
  checkmate::qassert(clusters, "S+")
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(observation, "S1")
  checkmate::qassert(value.y, "S1")
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(paired, "B1")
  checkmate::qassert(hide.ns, "B1")
  
  if (value.y != "percentage" && value.y != "absolute") {
    stop("The value.y name is invalid")
  }
  
  matrix.cell.count <- Celldata@matrix.cell.count
  
  if (!is.null(samples)) {
    matrix.cell.count <- matrix.cell.count[, colnames(matrix.cell.count) %in% samples]
  }
  
  metadata <- Celldata@metadata
  
  matrix.cell.count2 <- matrix.cell.count[rownames(matrix.cell.count) %in% clusters, ]
  cells.number <- sum(colSums(matrix.cell.count2))
  matrix.cell.count2 <- colSums(matrix.cell.count2)
  
  if (value.y == "percentage") {
    matrix.cell.count <- matrix.cell.count2 / base::apply(matrix.cell.count, 2, sum) * 100
  } else {
    matrix.cell.count <- matrix.cell.count2
  }
  
  matrix.cell.count <- reshape::melt(matrix.cell.count)
  
  matrix.cell.count <- base::merge(matrix.cell.count, metadata, by = "row.names")
  
  position_jitter <- ggplot2::position_jitter(seed = 42, width = 0.15)
  
  stat.test <- ggpubr::compare_means(data = matrix.cell.count,
                                     formula = stats::as.formula(paste0("value ~ ", observation)),
                                     method = test.statistics, paired = paired)
  stat.test <- data.frame(stat.test)
  
  stat.test$y.position <- max(matrix.cell.count$value)
  
  plot <- ggplot2::ggplot() +
    ggplot2::labs(title = "Boxplot representation",
                  subtitle = paste0("Clusters: ", paste0(clusters, collapse = ", "),
                                    " (", cells.number, " cells)")) +
    ggplot2::geom_boxplot(data = matrix.cell.count,
                          ggplot2::aes_string(x = observation,
                                              y = "value",
                                              fill = observation),
                          outlier.shape = NA, size = 0.2, fatten = 1) +
    ggiraph::geom_jitter_interactive(data = matrix.cell.count,
                                     ggplot2::aes_string(x = observation,
                                                         y = "value",
                                                         tooltip = "Row.names",
                                                         data_id = "individual"),
                                     color = "grey80", shape = 16, size = 1.5,
                                     position = position_jitter) +
    ggpubr::stat_pvalue_manual(stat.test, label = "p.signif", color = "#424242",
                               size = 5, hide.ns = hide.ns, step.increase = 0.1)
  
  if (value.y == "percentage") {
    plot <- plot +
      ggplot2::ylab("% abundance of cluster")
  } else {
    plot <- plot +
      ggplot2::ylab("absolute abundance of cluster")
  }
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "none")
  
  return(plot)
}

#' @title Plots of a volcano plot of statistical analysis
#'
#' @description This function aims to visualize the results of a differentially abundant a analysis using a Volcano plot.
#'
#' In such representation, each in dot corresponds to a cell cluster and dots are positioned in two dimensional space where the X-axis represents the log2(fold-change) and the Y-axis represents the -log10 of the p-value.
#' Un horizontal line is displayed accordingly to the p-value threshold and to vertical lines are displayed accordingly to the fold-change threshold.
#'
#' @param Celldata a Celldata object
#' @param comparison a character value containing the comparison to study
#' @param th.pv a numeric value containing the p-value threshold to use
#' @param th.fc a numeric value containing the fold-change threshold to use
#' @param plot.text a boolean value specifying if adds text directly at the plot
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotVolcano <- function(Celldata,
                        comparison,
                        th.pv = 1.3,
                        th.fc = 1.5,
                        plot.text = TRUE) {
  
  checkmate::qassert(comparison, "S1")
  checkmate::qassert(th.pv, "N1")
  checkmate::qassert(th.fc, "N1")
  
  stats <- Celldata@statistic
  stats <- stats[stats$comparison == comparison, ]
  
  stats$log10.pvalue <- -log10(stats$pvalue)
  max.fc <- max(abs(stats$lfc))
  
  stats$dir <- "ns"
  stats$dir[stats$log10.pvalue > th.pv & stats$lfc > log2(th.fc)] <- "up" #red
  stats$dir[stats$log10.pvalue > th.pv & stats$lfc < -log2(th.fc)] <- "down" #green
  
  clusters.diff <- stats[stats$log10.pvalue >= th.pv & abs(stats$lfc) > log(th.fc) / log(2), ]$clusters
  stats.diff <- stats[stats$clusters %in% clusters.diff, ]
  
  plot <- ggplot2::ggplot() +
    ggplot2::labs(title = "Volcano plot representation",
                  subtitle = comparison) +
    ggiraph::geom_point_interactive(data = stats,
                                    ggplot2::aes_string(x = "lfc",
                                                        y = "log10.pvalue",
                                                        fill = "dir",
                                                        tooltip = "clusters",
                                                        data_id = "clusters"),
                                    shape = 21, color = "black") +
    ggplot2::geom_hline(yintercept = th.pv, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(log2(th.fc), -log2(th.fc)), linetype = "dashed")
  
  if (plot.text == TRUE) {
    plot <- plot +
      ggrepel::geom_text_repel(data = stats.diff,
                               ggplot2::aes_string(x = "lfc",
                                                   y = "log10.pvalue",
                                                   label = "clusters"),
                               size = 3,
                               box.padding = grid::unit(0.35, "lines"),
                               point.padding = grid::unit(0.3, "lines"))
  }
  
  plot <- plot +
    ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(title = ""),
                               values = c("down" = "green", "up" = "red")) +
    ggplot2::scale_x_continuous(limits = c(-max.fc, max.fc), breaks = seq(-10, 10, 1)) +
    ggplot2::scale_y_continuous(breaks = c(seq(0, 10, 1), th.pv)) +
    ggplot2::xlab("log2(fold-change)") +
    ggplot2::ylab("-log10(p-value)")
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   legend.position = "bottom")
  
  return(plot)
}

#' @title Plots of a distogram of marker co-expression
#'
#' @description This function aims to visualize the pairwise co-expression between all markers using a distogram representation.
#' Each tile corresponds to the co-expression between two markers and is gradient-colored based on the Pearson or Spearman correlation
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifier of the cluster to use
#' @param method a character value indicating the name of correlation method to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotDistogram <- function(Celldata,
                          clusters = NULL,
                          method = c("pearson","spearman")) {
  
  checkmate::qassert(clusters, c("0", "S+"))
  
  matrix.exp <- Celldata@matrix.expression
  samples <- Celldata@samples
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(samples, cluster, matrix.exp)
  colnames(proj) <- c("samples", "clusters", colnames(matrix.exp))
  
  if (!is.null(clusters)) {
    proj <- proj[proj$clusters %in% clusters, ]
  }
  
  proj <- proj[, -c(1, 2)]
  proj <- stats::na.omit(proj)
  
  cormat <- round(stats::cor(proj, method = method), 2)
  dist <- stats::as.dist(1 - cormat)
  hc <- stats::hclust(dist)
  cormat <- cormat[hc$order, hc$order]
  cormat[upper.tri(cormat, diag = TRUE)] <- NA
  
  markers <- colnames(cormat)
  dimnames(cormat) <- NULL
  melted.cormat <- reshape2::melt(cormat)
  
  plot <- ggplot2::ggplot(data = melted.cormat,
                          ggplot2::aes_string(x = "Var1",
                                              y = "Var2",
                                              fill = "value")) +
    ggplot2::ggtitle("Distogram") +
    ggplot2::geom_tile(color = "white")
  
  plot <- plot +
    ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black",
                                  midpoint = 0, limit = c(-1, 1), na.value = "white",
                                  name = "Pearson correlation",
                                  oob = scales::oob_squish) +
    ggplot2::annotate(geom = "text",
                      x = seq(1, length(markers)),
                      y = seq(1, length(markers)),
                      angle = -45,
                      size = 4,
                      label = markers,
                      hjust = 1)
  
  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal",
      legend.key = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7,
                                                   barheight = 1,
                                                   title.position = "top",
                                                   title.hjust = 0.5))
  
  rownames(cormat) <- markers
  colnames(cormat) <- markers
  
  plot$cor <- cormat
  
  return(plot)
}

#' @title Plots of a scatter plot of marker co-expression
#'
#' @description This function aims to visualize co-expression between two markers using a scatter representation
#'
#' @param Celldata a Celldata object
#' @param marker1 a character value specifying the first marker to be visualized
#' @param marker2 a character value specifying the second marker to be visualized
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotScatter <- function(Celldata,
                        marker1,
                        marker2,
                        samples = NULL,
                        clusters = NULL) {
  
  checkmate::qassert(marker1, "S1")
  checkmate::qassert(marker2, "S1")
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(clusters, c("0", "S+"))
  
  matrix.exp <- Celldata@matrix.expression
  sample <- Celldata@samples
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(sample, cluster, matrix.exp)
  
  proj <- proj[, c(marker1, marker2, "sample", "cluster")]
  colnames(proj) <- c("marker1", "marker2", "samples", "clusters")
  
  if (!is.null(samples)) {
    proj <- proj[proj$samples %in% samples, ]
  }
  
  if (!is.null(clusters)) {
    proj <- proj[proj$clusters %in% clusters, ]
  }
  
  proj$value <- computeCellDensities(proj, n = 100)
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = proj,
                        ggplot2::aes_string(x = "marker1",
                                            y = "marker2",
                                            color = "value"))
  
  plot <- plot +
    ggplot2::scale_color_gradient(low = "yellow",
                                  high = "red",
                                  oob = scales::oob_squish) +
    ggplot2::xlab(marker1) +
    ggplot2::ylab(marker2)
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "none")
  
  return(plot)
}

#' @title Plots of phenotype of cell clusters using parallels coordinates
#'
#' @description This function aims to visualize the characteristics of cell clusters using parallels coordinates.
#' Each line in the representation corresponds to a biological sample for which marker/gene expressions are displayed.
#'
#' @param Celldata a Celldata object
#' @param condition.samples a character vector containing the variables to be studied for the samples. Possible values are: 'condition' or 'timepoint'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param clusters a character vector containing the identifiers of the clusters to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCoordinatesClusters <- function(Celldata,
                                    condition.samples = c("condition", "timepoint"),
                                    samples = NULL,
                                    clusters) {
  
  condition.samples <- match.arg(condition.samples)
  
  checkmate::qassert(condition.samples, "S1")
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(clusters, "S+")
  
  matrix.exp <- Celldata@matrix.expression
  sample <- Celldata@samples
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(sample, cluster, matrix.exp)
  colnames(proj) <- c("samples", "clusters", colnames(matrix.exp))
  
  melt.matrix <- reshape::melt(matrix.exp,variables="id")
  
  quantile <- plyr::ddply(melt.matrix, "variable",
                          function(x) {
                            quantile(x$value,
                                     probs = c(0.05, 0.95),
                                     na.rm = TRUE)})
  colnames(quantile) <- c("marker", "lower.bound", "upper.bound")
  
  if (!is.null(samples)) {
    proj <- proj[proj$samples %in% samples, ]
  }
  
  proj <- proj[proj$clusters %in% clusters, ]
  
  cells.number <- dim(proj)[1]
  
  proj <- plyr::ddply(proj, "samples", function(x) {
    x$samples <- NULL
    x$clusters <- NULL
    apply(x, 2, stats::median)
  })
  
  proj <- reshape::melt(proj, id.vars = "samples")
  
  means <- plyr::ddply(proj, "variable",
                       function(df) {
                         mean(df$value, na.rm = TRUE)})
  colnames(means) <- c("marker", "means")
  
  max.value <- max(c(proj$value, quantile$upper.bound), na.rm = TRUE)
  min.value <- min(c(proj$value, quantile$lower.bound), na.rm = TRUE)
  
  max.value <- max.value * (1 + sign(max.value) * 0.1)
  min.value <- min.value * (1 - sign(min.value) * 0.1)
  
  Celldata@metadata$samples <- rownames(Celldata@metadata)
  proj <- merge(proj, Celldata@metadata, by = "samples")
  
  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(paste0("Parallel coordinates - clusters: ",
                            paste0(clusters, collapse = ", "), " (", cells.number, " cells)")) +
    ggplot2::geom_ribbon(data = quantile,
                         ggplot2::aes_string(x = "marker",
                                             ymin = "lower.bound",
                                             ymax = "upper.bound"),
                         alpha = 0.1, group = 1, fill = "grey20") +
    ggiraph::geom_line_interactive(data = proj,
                                   ggplot2::aes_string(x = "variable",
                                                       y = "value",
                                                       group = "samples",
                                                       color = condition.samples,
                                                       tooltip = "samples",
                                                       data_id = "individual"),
                                   size = 0.5, alpha = 1) +
    ggplot2::geom_line(data = means,
                       ggplot2::aes_string(x = "marker",
                                           y = "means"),
                       group = 1, linetype = "dashed", size = 1)
  
  plot <- plot +
    ggplot2::xlab("markers") +
    ggplot2::ylab("expressions") +
    ggplot2::scale_y_continuous(limits = c(min.value, max.value),
                                breaks = round(seq(0, max.value, by = 1), 0))
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "bottom")
  
  return(plot)
}

#' @title Plots a LDA representation based cell cluster abundances
#'
#' @description This function aims to represent a Linear Discriminant Analysis representation based on cell cluster abundances
#'
#' @param Celldata a Celldata object
#' @param levels a character value containing the variable to be displayed. Possible values are: 'clusters' or 'samples'
#' @param ref.condition a character value providing the name of reference condition
#' @param condition a character value providing the name of the condition to be compared
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#'
#' @return xx
#'
#' @export
#'
plotLDA <- function(Celldata,
                    levels = c("predictions", "coefficients"),
                    ref.condition,
                    condition,
                    clusters = NULL) {
  
  checkmate::qassert(levels, "S1")
  checkmate::qassert(condition, "S+")
  checkmate::qassert(ref.condition, "S+")
  checkmate::qassert(clusters, c("0", "S*"))
  
  if (levels != "predictions" && levels != "coefficients") {
    stop("The levels name is invalid")
  }
  
  matrix.abundance <- Celldata@matrix.abundance
  
  if (!is.null(clusters)) {
    matrix.abundance <- matrix.abundance[rownames(matrix.abundance) %in%
                                           clusters, ]
  }
  
  values.condition <- matrix.abundance[, grepl(condition, colnames(matrix.abundance)),
                                       drop = TRUE]
  values.ref.condition <- matrix.abundance[, grepl(ref.condition, colnames(matrix.abundance)),
                                           drop = TRUE]
  
  grouping <- c(rep("condition", ncol(values.condition)),
                rep("ref.condition", ncol(values.ref.condition)))
  
  matrix.abundance <- cbind(values.condition, values.ref.condition)
  matrix.abundance <- t(matrix.abundance)
  
  matrix.abundance <- matrix.abundance[,apply(matrix.abundance,2,sum)!=0]
  
  model <- MASS::lda(x = matrix.abundance, grouping = grouping)
  predict <- stats::predict(model, newdata = matrix.abundance)
  
  if (levels == "predictions") {
    
    coeffs <- data.frame(predict$posterior)
    coeffs$sample <- rownames(coeffs)
    coeffs$class <- as.character(grouping)
    
    coeffs$predicted <- "predicted: ref.condition"
    coeffs$predicted[coeffs$condition > 0.5] <- "predicted: condition"
    
    coeffs$probability <- coeffs$condition
    coeffs$probability[coeffs$probability < 0.5] <- coeffs$ref.condition[coeffs$probability < 0.5]
    
    coeffs$predicted <- factor(coeffs$predicted, levels = c("predicted: ref.condition", "predicted: condition"))
    plot <- ggplot2::ggplot() +
      ggplot2::ggtitle("Linear discriminant analysis - class predictions") +
      ggplot2::geom_col(data = coeffs,
                        ggplot2::aes_string(x = "sample",
                                            y = "probability",
                                            fill = "class"),
                        color = "black")
    
    plot <- plot +
      ggplot2::xlab("") +
      ggplot2::ylab("classification probability") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
      ggplot2::scale_fill_manual(values = c("ref.condition" = "#87C7A9",
                                            "condition" = "#6890B9")) +
      ggplot2::facet_wrap(~predicted, nrow = 2, scale = "free_y")
    
    plot <- plot +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank())
    
  } else {
    
    coeffs <- data.frame(model$scaling)
    coeffs$clusters <- rownames(coeffs)
    coeffs$clusters <- factor(coeffs$clusters,
                              levels = coeffs$clusters[order(coeffs$LD1, decreasing = TRUE)])
    
    mean.condition <- apply(values.condition, 1, mean)
    mean.ref.condition <- apply(values.ref.condition, 1, mean)
    log2fc <- data.frame(log2fc = log2(mean.condition / mean.ref.condition))
    log2fc$clusters <- rownames(log2fc)
    
    coeffs <- merge(coeffs, log2fc, by = "clusters")
    
    plot <- ggplot2::ggplot() +
      ggplot2::ggtitle("Linear discriminant analysis - cluster coefficients") +
      ggplot2::geom_col(data = coeffs,
                        ggplot2::aes_string(x = "clusters",
                                            y = "LD1",
                                            fill = "log2fc"))
    
    plot <- plot +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
      ggplot2::scale_fill_gradientn(colours = c("green", "black", "red"),
                                    limits = c(-1, 1),
                                    oob = scales::oob_squish)
    
    plot <- plot +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    
  }
  
  return(plot)
}

#' @title Plot a marker for each cluster using parallel coordinates
#'
#' @description This function aims to visualise the expression of a given marker in each cluster using parallel coordinates.
#' Each line in the representation corresponds to a biological sample for which marker/gene expressions are displayed.
#'
#' @param Celldata a Celldata object
#' @param condition.samples a character vector containing the variables to be studied for the samples. Possible values are: 'condition' or 'timepoint'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param marlers a character vector containing the name of the markers to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCoordinatesMarkers <- function(Celldata,
                                   condition.samples = c("condition", "timepoint"),
                                   samples = NULL,
                                   clusters = NULL,
                                   markers) {
  
  condition.samples <- match.arg(condition.samples)
  
  checkmate::qassert(condition.samples, "S1")
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(clusters, c("0", "S+"))
  checkmate::qassert(markers, "S+")
  
  matrix.exp <- Celldata@matrix.expression
  sample <- Celldata@samples
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(sample, cluster, matrix.exp)
  colnames(proj) <- c("samples", "clusters", colnames(matrix.exp))
  
  melt.matrix <- reshape::melt(proj, id.vars = c("samples","clusters"))
  
  if (!is.null(samples)) {
    melt.matrix <- melt.matrix[melt.matrix$samples %in% samples, ]
  }
  if (!is.null(clusters)) {
    melt.matrix <- melt.matrix[melt.matrix$clusters %in% clusters, ]
  }
  
  melt.matrix <- melt.matrix[melt.matrix$variable %in% markers, ]
  
  proj2 <- plyr::ddply(melt.matrix, c("clusters","samples"), function(x) {
    stats::median(x$value)
  })
  colnames(proj2) <- c("clusters", "samples", "value")
  
  means <- plyr::ddply(proj2, "clusters",
                       function(df) {
                         mean(df$value, na.rm = TRUE)})
  colnames(means) <- c("clusters", "means")
  
  Celldata@metadata$samples <- rownames(Celldata@metadata)
  proj2 <- merge(proj2, Celldata@metadata, by = "samples")
  
  plot <- ggplot2::ggplot() +
    ggplot2::ggtitle(paste0("Parallel coordinates - Markers: ",
                            paste0(markers))) +
    ggiraph::geom_line_interactive(data = proj2,
                                   ggplot2::aes_string(x = "clusters",
                                                       y = "value",
                                                       group = "samples",
                                                       color = condition.samples,
                                                       tooltip = "samples",
                                                       data_id = "individual"),
                                   size = 0.5, alpha = 1)
  ggplot2::geom_line(data = means,
                     ggplot2::aes_string(x = "clusters",
                                         y = "means"),
                     group = 1, linetype = "dashed", size = 1)
  
  plot <- plot +
    ggplot2::xlab("clusters") +
    ggplot2::ylab("expressions")
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "bottom")
  
  return(plot)
}

#' @title Plot compare cluster 
#'
#' @description This function aims to visualise the expression of a given marker in each cluster.
#'
#' @param Celldata a Celldata object
#' @param clusters1 a character vector containing the identifier of the cluster to use
#' @param clusters2 a character vector containing the identifier of the cluster to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCompareClusters <- function(Celldata,
                                clusters1,
                                clusters2) {
  
  checkmate::qassert(clusters1, "S1")
  checkmate::qassert(clusters2, "S1")
  
  matrix.exp <- Celldata@matrix.expression
  cluster <- Celldata@identify.clusters
  
  proj <- cbind(matrix.exp, cluster)
  
  parameters <- colnames(proj)
  markers <- parameters[!parameters %in% c("cluster")]
  
  exp.values.clusters1 <- proj[proj$cluster %in% clusters1, ]
  exp.values.clusters2 <- proj[proj$cluster %in% clusters2, ]
  
  exp.values.clusters <- rbind(exp.values.clusters1, exp.values.clusters2)
  
  melt.matrix <- reshape::melt(exp.values.clusters, id.vars = c("cluster"))
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_density(data = melt.matrix,
                          mapping = ggplot2::aes_string(x = "value",
                                                        fill = "cluster"),
                          size = 0.5,
                          alpha = 0.8) +
    ggplot2::facet_wrap(~variable)
  
  plot <- plot +
    ggplot2::scale_fill_manual(values = c("#ff0000", "#0000ff")) +
    ggplot2::scale_x_continuous(limits = c(-2, 6), breaks = c(-2:6)) +
    ggplot2::scale_y_continuous(labels = function(x) {format(round(x, 1),nsmall = 1)})
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.35, "cm"))
  
  return(plot)
}
