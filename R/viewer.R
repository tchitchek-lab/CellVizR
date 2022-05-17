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
#' @param UMAPdata a UMAPdata object
#' @param stats a character vector providing the statistics to display. Possible values are: 'min', 'median', 'mean', 'q75', 'max'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param sort a boolean value indicating if clusters must be sorted by the number associated sample
#'  
#' @return a ggplot2 object
#' 
#' @export
#'    
plotCellCounts = function(UMAPdata, 
                          stats = c("min","median","mean","q75","max"),
                          samples = NULL,
                          sort = TRUE) {
  
  checkmate::qassert(stats, "S*")
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(sort, "B1")
  
  nb.cells = UMAPdata@matrix.expression
  nb.cells$samples = UMAPdata@samples
  
  if (!is.null(samples)) {
    nb.cells = nb.cells[nb.cells$samples %in% samples,]
  }
  
  nb.cells = plyr::ddply(nb.cells, "samples", nrow) 
  colnames(nb.cells) = c("samples", "nb")
  
  if (sort == TRUE) {
    nb.cells = nb.cells[order(nb.cells$nb, decreasing = TRUE),]
    nb.cells$samples = factor(nb.cells$samples, levels = nb.cells$samples)
  } else {
    nb.cells$samples <- factor(nb.cells$samples,levels= nb.cells$samples)
  }
  
  min = min(nb.cells$nb)
  mean = mean(nb.cells$nb)
  median = stats::median(nb.cells$nb)
  q75 = stats::quantile(nb.cells$nb, probs = 0.75)
  max = max(nb.cells$nb)
  
  xscale.middle = nrow(nb.cells)/2
  if(xscale.middle%%2==0){
    xscale.middle = xscale.middle+0.5
  }
  
  plot = ggplot2::ggplot() +
    ggplot2::geom_segment(data = nb.cells, ggplot2::aes_string(x="samples", xend="samples", yend="nb"), y=0, color = "gray50") + 
    ggplot2::geom_point(data = nb.cells, ggplot2::aes_string(x="samples", y="nb"), size=2)
  
  if("min" %in% stats){
    plot = plot +
      ggplot2::geom_hline(yintercept=min, color="blue") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="min"), label=paste0("min: ", format(min, scientific = FALSE, big.mark = ",")), color="blue", vjust=0)
  }
  
  if("max" %in% stats){
    plot = plot +
      ggplot2::geom_hline(yintercept=max, color="blue") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="max"), label=paste0("max: ", format(max, scientific = FALSE, big.mark = ",")), color="blue", vjust=0)
  }
  
  if("median" %in% stats){
    plot = plot +
      ggplot2::geom_hline(yintercept=median, color="red") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="median"),label=paste0("median: ", format(median, scientific = FALSE, big.mark = ",")), color="red", vjust=0)
  }
  
  if("mean" %in% stats){
    plot = plot +
      ggplot2::geom_hline(yintercept=mean, color="orange") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="mean"), label=paste0("mean: ", format(mean, scientific = FALSE, big.mark= "," )),color="orange", vjust=0)
  }
  
  if("q75" %in% stats){
    plot = plot +
      ggplot2::geom_hline(yintercept=q75, color="purple") +
      ggplot2::geom_text(data = nb.cells, ggplot2::aes_string(x="xscale.middle", y="q75"), label=paste0("q75: ", format(q75, scientific = FALSE, big.mark = ",")), color="purple", vjust=0)
  }
  
  plot = plot +
    ggplot2::scale_y_continuous(labels = function(x){format(x, scientific = FALSE, big.mark=",")}) +
    ggplot2::ylab("number of cells") 
  
  plot = plot + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"), 
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
                   axis.title.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), 
                   panel.border = ggplot2::element_rect(fill=NA), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(), 
                   legend.position = "none")
  
  invisible(plot) 
  
}

#' @title Plots the numbers of cells for each clusters
#'
#' @description This function aims to visualize the number of cells associated to each clusters.
#' 
#' This representation displays the clusters in the X-axis and the total number of associated cells in the Y-axis.
#' 
#' @param UMAPdata a UMAPdata object
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param sort a boolean value indicating if clusters must be sorted by the number associated cluster
#'  
#' @return a ggplot2 object
#' 
#' @export
#'    
plotClustersCounts = function(UMAPdata, 
                              clusters = NULL,
                              sort = TRUE) {
  
  checkmate::qassert(stats, "S*")
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(sort, "B1")
  
  matrix.cell = UMAPdata@matrix.cell.count
  
  if (!is.null(clusters)) {
    matrix.cell = matrix.cell[rownames(matrix.cell) %in% clusters,]
  } else {
    clusters = unique(UMAPdata@identify.clusters)
    matrix.cell = matrix.cell[clusters, , drop = FALSE]
  }
  
  cells.number = sum(colSums(matrix.cell))
  matrix.cell = cbind(matrix.cell, "sum.of.samples" = apply(matrix.cell, 1, sum))
  
  matrix.cell = cbind("clusters" = rownames(matrix.cell), matrix.cell)
  
  if (sort == TRUE) {
    order <- order(matrix.cell$sum.of.samples,decreasing=TRUE)
    matrix.cell$clusters <- factor(matrix.cell$clusters,levels=matrix.cell$clusters[order])
  } else {
    matrix.cell$clusters <- factor(matrix.cell$clusters,levels=clusters)
  }
  
  data.melted <- reshape2::melt(matrix.cell, id = c("clusters"))
  colnames(data.melted) <- c("clusters", "samples", "values")
  data.melted$total <- ifelse(data.melted[, "samples"] == "sum.of.samples", "sum for selected samples", "")
  
  plot = ggplot2::ggplot() + 
    ggplot2::ggtitle(paste0("Count viewer (", format(cells.number, big.mark = " "), " cells)", sep="")) + 
    ggplot2::geom_jitter(data = subset(data.melted, samples != "sum.of.samples"),
                         ggplot2::aes_string(x="clusters", y="values", size="values", fill="samples"),
                         height = 0, width = 0.5, shape = 21, alpha = 0.4) + 
    ggplot2::geom_point(data = subset(data.melted, samples == "sum.of.samples"),
                        ggplot2::aes_string(x="clusters", y="values", size="values", shape = "total"),
                        fill = "grey40")
  
  plot = plot + 
    ggplot2::scale_size_area(name = "Number of cells") + 
    ggplot2::scale_shape_manual(values = 21) + 
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 1.1 * max(data.melted$values))) + 
    ggplot2::xlab("clusters") + 
    ggplot2::ylab("Number of cells")
  
  plot = plot + 
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), 
      legend.key = ggplot2::element_blank(), 
      plot.title = ggplot2::element_text(hjust = 0.5))
  
  invisible(plot)
  
}

#' @title Plots of phenotype of identified cell clusters 
#'
#' @description This function aims to visualize xxx
#' 
#' @param UMAPdata a UMAPdata object
#' @param cluster a character vector containing the identifier of the cluster to use
#' @param quant.low a numeric value containg xxx
#' @param quant.high a numeric value containing xxx
#' @param dip.th a numeric value containing xxx
#'  
#' @return a ggplot2 object
#' 
#' @export
#'    
plotPhenoClusters = function(UMAPdata,
                             cluster,
                             quant.low = 0.05,
                             quant.high = 0.95, 
                             dip.th = 0.01) {
  
  checkmate::qassert(cluster, "N+")
  checkmate::qassert(quant.low, "N1")
  checkmate::qassert(quant.high, "N1")
  checkmate::qassert(dip.th, "N1")
  
  expression_color_palette = c("white","yellow","orange","red","red4")
  
  matrix.exp = UMAPdata@matrix.expression
  samples = UMAPdata@samples
  clusters = UMAPdata@identify.clusters
  assignments = UMAPdata@metadata
  
  proj = cbind(samples, matrix.exp, clusters)
  
  parameters = colnames(proj)
  markers = parameters[!parameters %in% c("samples", "clusters")]
  
  grob = list()
  
  for (marker in markers) {
    
    exp.values = proj[,c(marker, "clusters", "samples")]
    colnames(exp.values) = c("exp", "clusters", "samples")
    quantiles = stats::quantile(exp.values$exp, probs = c(quant.low, quant.high))
    
    print(cluster)
    exp.values.clusters = exp.values[exp.values$clusters %in% cluster,]
    
    order = unique(assignments$timepoint)
    assignments = assignments[samples, , drop = FALSE]
    exp.values$timepoint = assignments[exp.values$samples,"timepoint"]
    order = intersect(order, unique(assignments$timepoint))
    exp.values$timepoint = factor(exp.values$timepoint, levels = order)
    
    median = median(exp.values.clusters$exp)
    tmp.exp = exp.values.clusters$exp
    dip = diptest::dip.test(tmp.exp)
    
    if(dip$p.value<dip.th){
      color = "red"
    }else{
      color = "green"
    }
    
    plot = ggplot2::ggplot() + 
      ggplot2::ggtitle(marker) + 
      ggridges::geom_density_ridges_gradient(data = exp.values, ggplot2::aes(x=exp, fill=stat(x), y=0)) + 
      ggplot2::scale_fill_gradientn(limits=c(quantiles[1], quantiles[2]),
                                    colours = expression_color_palette, na.value = "black") +
      ggplot2::geom_density(data = exp.values.clusters, ggplot2::aes_string(x="exp"), 
                            color = color, fill=NA, size=0.75) + 
      # ggplot2::geom_density(data = exp.values, ggplot2::aes_string(x="exp", color = "TP"), 
      #                       size = 0.75) + group = "samples",
      ggplot2::geom_vline(xintercept = median, linetype = "dashed", color = "gray20") +
      ggplot2::scale_x_continuous(limits = c(-2,6), breaks = c(-2:6)) +
      ggplot2::scale_y_continuous(labels = function(x) {format(round(x, 1), nsmall = 1)}) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        legend.key.width = grid::unit(0.35, "cm"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    grob[[marker]] = plot
  }
  
  test = gridExtra::grid.arrange(grobs = grob,
                                 ncol = 6)
  
  return(test)
  
}

#' @title Plots a representation of a computed manifold
#'
#' @description This function aims to visualize a computed manifold representation for given analysis.
#' 
#' This representation can be used on a UMAPdata object for which a manifold analysis has been performed.
#' 
#' If a cell clustering has been performed, then the clusters are delineated using concave hulls.
#' Additionnaly, the manifold can be colored based on the local cell density or marker expressions. 
#' It is possible to centred ans reduce the values of expressions. 
#'  
#' @param UMAPdata a UMAPdata object
#' @param markers a character value providing the name of the marker to use for the colouring. By default, cells are colored based on their local density
#' @param samples a character vector containing the names of biological samples to use
#' @param scale a boolean value specifing if expression calue must be rescaled
#' @param quant.low a numeric value providing the number of first quantile
#' @param quant.high a numeric value providing the number of last quantile  
#'  
#' @return a ggplot2 object
#' 
#' @export
#'    
plotManifold = function(UMAPdata, 
                        markers = "density",
                        samples = NULL,
                        scale = FALSE,
                        quant.low = 0.05,
                        quant.high = 0.95) {
  
  checkmate::qassert(markers, "S1")
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(quant.low, "N1")
  checkmate::qassert(quant.high, "N1")
  
  proj = UMAPdata@manifold
  colnames(proj) = c("UMAP1","UMAP2")
  clusters = UMAPdata@identify.clusters
  
  if(length(clusters)!=0){
    proj$clusters = clusters
  }
  
  proj = na.omit(proj)
  
  lines = UMAPdata@concave.hulls
  
  if(markers=="density"){
    legend.title =  "cell density"
    proj$value = computeCellDensities(proj, n=100)
  }else{
    legend.title =  paste0(markers," expression")
    proj$value = UMAPdata@matrix.expression[,markers]
    palette = rev(grDevices::rainbow(100,alpha=1)[seq(1,85)])
  }
  
  plot = ggplot2::ggplot()
  
  if(!is.null(samples)){
    proj.ref = proj
    proj     = proj[UMAPdata@samples %in% samples,]
    plot = plot +
      ggplot2::geom_point(data=proj.ref, ggplot2::aes_string(x="UMAP1", y="UMAP2"), color="gray", size=0.0001) +
      ggnewscale::new_scale_color()
  }
  
  plot = plot +
    ggplot2::geom_point(data=proj, ggplot2::aes_string(x="UMAP1", y="UMAP2", color="value"), size=0.0001)
  
  
  if(nrow(lines)!=0){
    plot = plot + ggplot2::geom_path(data=lines, ggplot2::aes_string(x="dim1", y="dim2", group="clusters"), color="black", size=0.2)
  }
  
  if(markers=="density"){
    plot = plot + ggplot2::scale_color_gradientn(colours=c("black","red"),
                                                 labels=function(x){return(rep("", length(x)))})
  }else{
    limits = abs.proj(proj, scale, quant.low, quant.high)
    plot = plot + ggplot2::scale_color_gradientn(colours = palette,limits = limits)
  }
  
  plot = plot+
    ggplot2::labs(title = "UMAP representation", color=legend.title)
  
  plot = plot+
    ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01,0.01))
  
  plot = plot + 
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   aspect.ratio = 1,
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(hjust = 1, vjust = 1, face = 'bold'))
  
  return(plot)
  
}

#' @title Plots a PCA representation based cell cluster abundances 
#'
#' @description This function aims to represent a Principal Component Analysis representation based on cell cluster abundances.
#' In such representation, clusters or samples are positioned based on computed principal components. 
#' The representation can be displayed based on specific principal components. 
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples. 
#' 
#' @param UMAPdata a UMAPdata object
#' @param levels a character value containing the variable to be displayed. Possible values are: 'clusters' or 'samples'
#' @param clusters a character vector containing the identifier of the cluster to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param components a numeric vector providing the components to display 
#' @param cor.radius.th a numeric value containing
#' @param cluster.label.size a numeric value containing xxx
#' @param sample.label.size a numeric value containing xxx
#' @param cluster.dot.size a numeric value containing xxx
#' @param sample.dot.size a numeric value containing xxx
#'     
#' @return a ggplot2 object
#' 
#' @export
#' 
plotPCA <- function(UMAPdata,
                    levels = c("clusters", "samples", "both"),
                    clusters = NULL,
                    samples = NULL,
                    components = c(1,2),
                    cor.radius.th = 0.6,
                    cluster.label.size = 3,
                    sample.label.size = 3,
                    cluster.dot.size = 2,
                    sample.dot.size = 3){
  
  levels <- match.arg(levels)
  
  checkmate::qassert(levels, "S1")
  checkmate::qassert(clusters, c("0","S*"))
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(components, "N2")
  checkmate::qassert(cor.radius.th, "N1")
  checkmate::qassert(cluster.label.size, "N1")
  checkmate::qassert(sample.label.size, "N1")
  checkmate::qassert(cluster.dot.size, "N1")
  checkmate::qassert(sample.dot.size, "N1")
  
  circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
    r <- diameter / 2
    t <- seq(0,2*pi,length.out = npoints)
    x <- center[1] + r * cos(t)
    y <- center[2] + r * sin(t)
    return(data.frame(x = x, y = y))
  }
  
  matrix.abundance = UMAPdata@matrix.abundance
  
  if(!is.null(clusters)) {
    matrix.abundance = matrix.abundance[rownames(matrix.abundance) %in% clusters, ]
  }
  if(!is.null(samples)) {
    matrix.abundance = matrix.abundance[, colnames(matrix.abundance) %in% samples]
  }
  
  res.PCA         = FactoMineR::PCA(t(matrix.abundance), graph=FALSE)
  var.explained   = res.PCA$eig[,2]
  
  data_clusters    = data.frame(res.PCA$var$coord[,components])
  colnames(data_clusters) = c("x","y")
  data_clusters$id = rownames(matrix.abundance)
  
  data_variables = data.frame(res.PCA$ind$coord[,components])
  colnames(data_variables) = c("x","y")
  data_variables$id = colnames(matrix.abundance)
  
  data_variables$x = scales::rescale(data_variables$x,to=c(-1,1))
  data_variables$y = scales::rescale(data_variables$y,to=c(-1,1))
  
  circle1 <- circleFun(c(0,0),2,npoints = 100)
  circle2 <- circleFun(c(0,0),2*cor.radius.th,npoints = 100)
  
  data_clusters = data_clusters[sqrt(data_clusters$x^2+data_clusters$y^2)>cor.radius.th,]
  
  xlab = paste0("PCA",components[1]," (",round(var.explained[components[1]],2),"%)")
  ylab = paste0("PCA",components[2]," (",round(var.explained[components[2]],2),"%)")
  
  plot <- ggplot2::ggplot()
  
  if (levels == "clusters" || levels == "both") {
    
    plot <- plot+
      ggplot2::geom_hline(yintercept=0,linetype = "dashed",size=0.2) + 
      ggplot2::geom_vline(xintercept=0,linetype = "dashed",size=0.2) + 
      ggplot2::geom_path(data=circle1,ggplot2::aes_string(x="x",y="y"),size=0.2,color="gray") +
      ggplot2::geom_path(data=circle2,ggplot2::aes_string(x="x",y="y"),size=0.2,linetype="dashed",color="gray") +
      ggplot2::geom_point(data = data_clusters, ggplot2::aes_string(x="x", y="y"),shape=21,size=cluster.dot.size,stroke=0.1,fill="black",color="black") +
      ggrepel::geom_text_repel(data = data_clusters, ggplot2::aes_string(x="x", y="y", label="id"),
                               size = cluster.label.size,
                               color = "black",
                               max.overlaps = Inf,
                               min.segment.length = 0,
                               segment.color = NA,
                               segment.size = 0.1)
    
  }
  
  if (levels == "samples" || levels == "both") {
    
    plot <- plot +
      ggplot2::geom_hline(yintercept=0,linetype = "dashed",size=0.2) + 
      ggplot2::geom_vline(xintercept=0,linetype = "dashed",size=0.2) + 
      ggplot2::geom_point(data = data_variables, ggplot2::aes_string(x="x", y="y"),shape=21,size=sample.dot.size,stroke=0.1,fill="gray",color="black") +
      ggrepel::geom_text_repel(data = data_variables, ggplot2::aes_string(x="x", y="y", label="id"),
                               size=sample.label.size,
                               color = "black",
                               max.overlaps = Inf,
                               min.segment.length = 0,
                               segment.color = NA,
                               segment.size = 0.1)
  }
  
  plot <- plot + ggplot2::labs(title="PCA represnetation")
  
  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(color=NA),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(hjust = 0.5))
  
  
  plot
  
  
  return(plot)
}

#' @title Plots a MDS representation based on cell cluster abundances
#'
#' @description This function aims to visualize the similarities between samples or clusters based on their abundances, using a Multidimensional Scaling representation.
#' Each dot represents a sample or a cluster and the distances between the dots are proportional to the Euclidean distance between these objects.
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples. 
#'     
#' @param UMAPdata a UMAPdata object
#' @param levels a character value containing the variable to be displayed. Possible values are: 'clusters' or 'samples'
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used 
#'  
#' @return a ggplot2 object
#' 
#' @export
#' 
plotMDS <- function(UMAPdata,
                    levels = c("clusters", "samples"),
                    clusters = NULL,
                    samples = NULL){
  
  levels <- match.arg(levels)
  
  checkmate::qassert(levels, "S1")
  checkmate::qassert(clusters, c("0","S*"))
  checkmate::qassert(samples, c("0","S*"))

  matrix.abundance = UMAPdata@matrix.abundance
  
  if (levels == "clusters") { 
    
    if(!is.null(clusters)) {
      matrix.abundance = matrix.abundance[rownames(matrix.abundance) %in% clusters, ]
    }
    
    dist.clusters <- stats::dist(matrix.abundance)
    
    fit1 <- MASS::isoMDS(dist.clusters, k=2) 
    x1 <- fit1$points[,1]
    y1 <- fit1$points[,2]
    
    min.lim = min(min(x1),min(y1))*1.1
    max.lim = max(max(x1),max(y1))*1.1
    
    proj.clusters = data.frame(x=x1,
                               y=y1)
    
    proj.clusters$clusters = rownames(matrix.abundance)
    
    plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ",round(fit1$stress,2))) +
      ggplot2::geom_hline(yintercept=(min.lim+max.lim)/2, linetype = "dashed") + 
      ggplot2::geom_vline(xintercept=(min.lim+max.lim)/2, linetype = "dashed") + 
      #geom_polygon(data = hulls, aes(x=x,y=y,group=cond,fill=color), alpha = 0.6) +
      ggplot2::geom_point(data = proj.clusters, ggplot2::aes_string(x="x", y="y",fill="clusters"),shape=21,color="black") +
      ggrepel::geom_text_repel(data = proj.clusters, ggplot2::aes_string(x="x", y="y",label="clusters"),size=4) 
    
    plot <- plot + ggplot2::coord_fixed()
    
    plot <- plot + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5,face="italic"),
                     aspect.ratio = 1,
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(), 
                     axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill=NA),
                     panel.grid.minor=ggplot2::element_blank(),
                     panel.grid.major=ggplot2::element_blank(),
                     legend.position="none")
    
  } else if (levels == "samples") { 
    
    if (!is.null(samples)) {
      matrix.abundance = matrix.abundance[, names(matrix.abundance) %in% samples]
    } 
    
    dist.samples <- stats::dist(t(matrix.abundance))
    
    fit2 <- MASS::isoMDS(dist.samples, k=2) 
    x2 <- fit2$points[,1]
    y2 <- fit2$points[,2]
    
    min.lim = min(min(x2),min(y2))*1.1
    max.lim = max(max(x2),max(y2))*1.1
    
    proj.samples = data.frame(x=x2,
                              y=y2)
    
    proj.samples$samples = colnames(matrix.abundance)
    
    plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ",round(fit2$stress,2))) +
      ggplot2::geom_hline(yintercept=(min.lim+max.lim)/2, linetype = "dashed") + 
      ggplot2::geom_vline(xintercept=(min.lim+max.lim)/2, linetype = "dashed") + 
      #geom_polygon(data = hulls, aes(x=x,y=y,group=cond,fill=color), alpha = 0.6) +
      ggplot2::geom_point(data = proj.samples, ggplot2::aes_string(x="x", y="y",fill="samples"),shape=21,color="black") +
      ggrepel::geom_text_repel(data = proj.samples, ggplot2::aes_string(x="x", y="y",label="samples"),size=4) 
    
    plot <- plot + ggplot2::coord_fixed()
    
    plot <- plot + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5,face="italic"),
                     aspect.ratio = 1,
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(), 
                     axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill=NA),
                     panel.grid.minor=ggplot2::element_blank(),
                     panel.grid.major=ggplot2::element_blank(),
                     legend.position="none")
  }
  
  return(plot)
}

#' @title Plots cell cluster abundances using a boxplot representation 
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#' 
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed. 
#' The representation can be resticted to a specific set of samples. 
#' Moreover boxplot can be constructed based on sample meta information. 
#' Statistic can be computed for all comparisons. 
#' 
#' @param UMAPdata a UMAPdata object
#' @param clusters a character vector containing the identifiers of the clusters to use
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observation a character value containing the parameters to use
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied 
#' @param hide.ns a boolean value indicating xxx
#'  
#' @return a ggplot2 object
#' 
#' @export
#' 
plotBoxplot = function(UMAPdata,
                       clusters, 
                       samples = NULL, 
                       observation = c("individual", "condition", "timepoint"),
                       test.statistics = c("wilcox.test", "t.test"),
                       paired = FALSE,
                       hide.ns = TRUE){
  
  observation <- match.arg(observation)
  test.statistics <- match.arg(test.statistics)
  
  checkmate::qassert(clusters, "N+")
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(observation, "S1")
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(paired, "B1")
  checkmate::qassert(hide.ns, "B1")

  matrix.cell.count = UMAPdata@matrix.cell.count
  
  if(!is.null(samples)) {
    matrix.cell.count = matrix.cell.count[, colnames(matrix.cell.count) %in% samples]
  }
  
  metadata = UMAPdata@metadata
  
  matrix.cell.count = matrix.cell.count[rownames(matrix.cell.count) %in% clusters,] 
  matrix.cell.count = base::apply(matrix.cell.count, 2, sum)
  matrix.cell.count = matrix.cell.count/base::apply(UMAPdata@matrix.cell.count, 2, sum)*100
  
  matrix.cell.count = reshape::melt(matrix.cell.count)
  
  matrix.cell.count = merge(matrix.cell.count, metadata, by = "row.names")
  
  position_jitter = ggplot2::position_jitter(seed=42,width=0.15)
  
  stat.test  = ggpubr::compare_means(data=matrix.cell.count, formula= as.formula(paste0("value ~ ",observation)), method = test.statistics, paired = paired)
  stat.test  = data.frame(stat.test)
  stat.test$y.position = max(matrix.cell.count$value)
  
  
  plot = ggplot2::ggplot() +
    ggplot2::ggtitle(paste0("abundance of cluster: ", paste0(clusters, collapse = ", "))) +
    ggplot2::geom_boxplot(data=matrix.cell.count,ggplot2::aes_string(x=observation,y="value"),
                          fill="gray",outlier.shape=NA,size=0.2,fatten = 1) +
    ggplot2::geom_jitter(data=matrix.cell.count,ggplot2::aes_string(x=observation,y="value"),
                         fill="black",shape=21,size=0.1,position=position_jitter) +
    ggpubr::stat_pvalue_manual(stat.test,label = "p.signif",color = "purple",size=5,hide.ns=hide.ns)
  
  
  plot = plot+
    ggplot2::ylab("abundance of cluster relative to parent gate")
  
  plot = plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust=0.5),
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
#'     
#' @param UMAPdata a UMAPdata object
#' @param th.pv a numeric value containing the p-value threshold to use
#' @param th.fc a numeric value containing the fold-change threshold to use 
#' 
#'  
#' @return a ggplot2 object
#' 
#' @export
#' 
plotVolcanoPlot = function(UMAPdata,
                           th.pv = 1.3,
                           th.fc = 1.5){
  
  checkmate::qassert(th.pv, "N1")
  checkmate::qassert(th.fc, "N1")
  
  stats = UMAPdata@statistic
  stats$log10.pvalue = -log10(stats$pvalue)
  max.fc = max(abs(stats$lfc))
  
  stats$dir = "none"
  stats$dir[stats$log10.pvalue > th.pv & stats$lfc > th.fc] = "red"
  stats$dir[stats$log10.pvalue > th.pv & stats$lfc < th.fc] = "green"
  
  plot = ggplot2::ggplot() +
    ggplot2::geom_point(data=stats, ggplot2::aes_string(x="lfc",y="log10.pvalue",fill="dir"),
                        shape=21,color="black") +
    ggplot2::scale_fill_manual(values=c("none"="gray","green"="green","red"="red")) +
    ggplot2::scale_x_continuous(limits=c(-max.fc,max.fc),breaks=seq(-10,10,1)) +
    ggplot2::scale_y_continuous(breaks=c(seq(0,10,1), th.pv)) +
    ggplot2::geom_hline(yintercept= th.pv, linetype = "dashed") + 
    ggplot2::geom_vline(xintercept=log2(th.fc), linetype = "dashed") + 
    ggplot2::geom_vline(xintercept=-log2(th.fc), linetype = "dashed") + 
    ggplot2::theme_bw() +
    ggplot2::xlab("log2(fold-change)") + 
    ggplot2::ylab("-log10(p-value)") +
    ggplot2::theme(legend.position = "none")
  
  return(plot) 
  
  
}

#' @title Plots of marker co-expression
#'
#' @description This function aims to visualize xxx
#' 
#' @param UMAPdata a UMAPdata object
#' @param clusters a character vector containing the identifier of the cluster to use
#' @param show.on.device a numeric value containg xxx
#'  
#' @return a ggplot2 object
#' 
#' @export
#'    
plotDistogram = function(UMAPdata,
                         clusters = NULL,
                         show.on.device = TRUE) {
  
  
  matrix.exp = UMAPdata@matrix.expression
  samples = UMAPdata@samples
  cluster = UMAPdata@identify.clusters
  
  proj = cbind(samples, cluster, matrix.exp)
  
  if(is.null(clusters)) {
    clusters = unique(cluster)
  } else if (all(clusters %in% UMAPdata@identify.clusters)) {
    if (typeof(clusters) != "character") {
      stop("Error in plotDistogram : 'clusters' parameter must be character vector")
    }
    
    clusters = unique(clusters)
    clusters.select = proj[, "cluster"] %in% clusters
    proj = proj[clusters.select,]
    
  } else {
    stop("Error in plotDistogram:\nUnknown clusters : ", paste0(setdiff(unique(clusters), unique(cluster)), collapse = " "))
  }
  
  proj = proj[, -c(1,2)]
  proj = stats::na.omit(proj)
  
  cormat = round(stats::cor(proj, method = "pearson"), 2)
  dist = stats::as.dist(1 - cormat)
  hc = stats::hclust(dist)
  cormat = cormat[hc$order, hc$order]
  cormat[upper.tri(cormat, diag = TRUE)] <- NA
  
  markers = colnames(cormat)
  dimnames(cormat) = NULL
  melted.cormat = reshape2::melt(cormat)
  
  plot = ggplot2::ggplot(data = melted.cormat, ggplot2::aes_string(x="Var1", y="Var2", fill="value")) +
    ggplot2::ggtitle("Distogram") + 
    ggplot2::geom_tile(color = "white") + 
    ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black",
                                  midpoint = 0, limit = c(-1,1), na.value = "white",
                                  name = "Pearson correlation") + 
    ggplot2::annotate(geom = "text", 
                      x = 1:length(markers),
                      y = 1:length(markers),
                      angle = -45,
                      size = 4,
                      label = markers,
                      hjust = 1) + 
    ggplot2::coord_fixed() +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   legend.justification = c(1, 0),
                   legend.position = c(0.4, 0.7),
                   legend.direction = "horizontal",
                   legend.key = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust=0.5)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7,
                                                   barheight      = 1,
                                                   title.position = "top",
                                                   title.hjust    = 0.5))
  
  if (show.on.device) {
    plot(plot)
  }
  
  rownames(cormat) = markers
  colnames(cormat) = markers
  
  plot$cor = cormat
  
  invisible(plot)
  
}

#' @title Plots of scatter
#'
#' @description This function aims to visualize xxx
#' 
#' @param UMAPdata a UMAPdata object
#' @param marker1 a character value containing xxx
#' @param marker2 a character value containing xxx
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param clusters a character vector containing the identifiers of the clusters to use. By default, all clusters are used
#'  
#' @return xx
#' 
#' @export
#'    
plotScatter = function(UMAPdata,
                       marker1,
                       marker2,
                       samples = NULL,
                       clusters = NULL) {
  
  matrix.exp   = UMAPdata@matrix.expression
  samples.xxx  = UMAPdata@samples
  clusters.xxx = UMAPdata@identify.clusters
  
  proj = cbind(samples.xxx, clusters.xxx, matrix.exp)
  
  proj = proj[,c(marker1,marker2,"samples.xxx","clusters.xxx")]
  colnames(proj) = c("marker1","marker2","samples","clusters")
  
  if(!is.null(samples)){
    proj = proj[proj$samples %in% samples]
  }
  if(!is.null(clusters)){
    proj = proj[proj$clusters %in% clusters]
  }
  
  proj$value = computeCellDensities(proj, n=100)
  
  plot = ggplot2::ggplot() +
    ggplot2::geom_point(data=proj,ggplot2::aes_string(x="marker1",y="marker2",color="value"))
  
  
  plot = plot+
    ggplot2::scale_color_gradient(low="yellow",high="red") + 
    ggplot2::xlab(marker1)+
    ggplot2::ylab(marker2)
  
  plot = plot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust=0.5),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "none")
  
  
  plot
  
}

#' @title Plots parallel coordinates
#'
#' @description This function aims to visualize xxx
#' 
#' @param UMAPdata a UMAPdata object
#'  
#' @return xx
#' 
#' @export
#'    
plotCoordinates = function(UMAPdata) {
  
}

#' @title Plots a LDA representation based cell cluster abundances 
#'
#' @description This function aims to represent a Linear Discriminant Analysis representation based on cell cluster abundances.
#' 
#' @param UMAPdata a UMAPdata object
#'  
#' @return xx
#' 
#' @export
#'    
plotLDA = function(UMAPdata) {
  
}

