#' @title Identify cell cluster of having similar marker expressions
#'
#' @description This function aims to identify cell clusters, which are groups of cells having similar expressions for selected markers, using different unsupervised clustering methods.
#'
#' Several clustering methods are available such as kmeans, kmedian, clara, DBSCAN, HDBSCAN and FlowSOM.
#' The cell clustering can be performed on the manifold representation or based on marker expression.
#'
#' @details
#' For each identify cell cluster, the boundaries of cells belonging to this cluster are delineated using a concave hull
#'
#' @param Celldata a Celldata object
#' @param space a character value containing the space of the clustering method to use. Possible values are: 'manifold' or 'markers'
#' @param markers a character vector providing the cell markers to use for the manifold generation
#' @param method a character value containing the name of the clustering method to use. Possible values are: 'kmeans', 'kmedian', 'clara', 'DBSCAN' and 'FlowSOM'
#' @param concavity a numeric value providing a relative measure of concavity for the computation of the concave hulls (please refer to the function 'concaveman' of the 'concaveman' package)
#' @param length.threshold a numeric value providing a threshold of the segment length for the computation of the concave hulls (please refer to the function 'concaveman' of the concaveman package)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param ... other arguments passed on to the methods
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#'
identifyClusters <- function(Celldata,
                             space = c("manifold", "markers"),
                             markers = NULL,
                             method = c("kmeans", "kmedian", "clara", "DBSCAN", "FlowSOM"),
                             concavity = 2,
                             length.threshold = 0,
                             seed = 42,
                             ...) {
  
  space <- match.arg(space)
  method <- match.arg(method)
  
  checkmate::qassert(space, "S1")
  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(method, "S1")
  checkmate::qassert(concavity, "N1")
  checkmate::qassert(length.threshold, "N1")
  checkmate::qassert(seed, "N1")
  
  proj  <- Celldata@manifold
  exprs <- Celldata@matrix.expression
  
  message(paste0("Clustering method is: ", method))
  message()
  message("Identifying cell clusters...")
  
  if (space == "manifold") {
    data <- proj
    clustering.markers <- Celldata@manifold.params$markers
  } else if (space == "markers") {
    if (is.null(markers)) {
      markers <- names(exprs)
    }
    data <- exprs[, names(exprs)
                  %in% markers]
    clustering.markers <- markers
  }
  
  do.call("set.seed", list(seed))
  switch(method,
         kmeans = {
           clusters <- stats::kmeans(data, ...)$cluster
         },
         kmedian = {
           clusters <- Gmedian::kGmedian(data, ...)$cluster[, 1]
         },
         clara = {
           clusters <- cluster::clara(data, ...)$clustering
         },
         DBSCAN = {
           clusters <- dbscan::dbscan(data, ...)$cluster
         },
         FlowSOM = {
           Flow <- FlowSOM::FlowSOM(as.matrix(data), colsToUse = clustering.markers, ...)
           clusters <- FlowSOM::GetMetaclusters(Flow)
         })
  
  clusters <- as.vector(clusters)
  clusters <- as.character(clusters)
  
  Celldata@identify.clusters <- clusters
  
  identify.Clusters.params <- list(...,
                                   concavity = concavity,
                                   length.threshold = length.threshold,
                                   clustering.markers = clustering.markers)
  Celldata@identify.clusters.params <- identify.Clusters.params
  
  if (space == "manifold") {
    message("computing cell clusters boundaries...")
    Celldata@concave.hulls <- computeConcaveHulls(proj = proj,
                                                  clusters = clusters,
                                                  concavity = concavity,
                                                  length.threshold = length.threshold)
  }
  
  samples <- Celldata@samples
  
  message("computing cell cluster count matrix...")
  Celldata@matrix.cell.count <- computeCellCounts(proj = data,
                                                  samples = samples,
                                                  clusters = clusters)
  
  # Renaming samples in Celldata@samples slot to accomodate the modifications made by computeCellCounts
  
  for(s in colnames(Celldata@matrix.cell.count))
  {
    replaceValueID = which(colnames(Celldata@matrix.cell.count) == s)
    
    replaceValue = unique(colnames(Celldata@matrix.cell.count))[replaceValueID]
    originalValue = unique(Celldata@samples)[replaceValueID]
    
    Celldata@samples[Celldata@samples == originalValue] = replaceValue
    
    
  }
  
  
  message("computing cell cluster abundance matrix...")
  count <- Celldata@matrix.cell.count
  Celldata@matrix.abundance <- computeClusterAbundances(count = count)
  
  
  # Renaming samples in Celldata@matrix.abundance slot to accomodate the modifications made by computeClusterAbundances
  
  colnames(Celldata@matrix.abundance) = colnames(Celldata@matrix.cell.count)
  
  # Renaming samples in Celldata@metadata slot to accomodate the modifications made by the previous tasks
  
  rownames(Celldata@metadata) = colnames(Celldata@matrix.cell.count)
  
  validObject(Celldata)
  
  return(Celldata)
}

# @title Internal - Computes the concave hulls for identified cell clusters
#
# @description This function is used internally to computes the concave hulls each cell cluster.
# Each concave hulls correspond the boundaries of the cells cluster with manifold representation.
#
# @param proj a data.frame providing the manifold representation
# @param clusters a character vector providing id of cell clusters for which the concave hulls must be computed
# @param concavity a numeric value providing a relative measure of concavity (please refer to the function 'concaveman' of the 'concaveman' package)
# @param length.threshold a numeric value providing a threshold for the segment length (please refer to the function 'concaveman' of the concaveman package)
#
# @return a data.frame containing the concave hulls for each cluster (dim1, dim2, clusters)
#
computeConcaveHulls <- function(proj,
                                clusters,
                                concavity = 2,
                                length.threshold = 0) {
  
  proj <- cbind(proj, clusters)
  colnames(proj) <- c("dim1", "dim2", "clusters")
  
  concave.hulls <- data.frame()
  for (k in sort(unique(proj$clusters))) {
    sub.exprs <- proj[proj$clusters == k, ]
    xcoord <- sub.exprs$dim1
    ycoord <- sub.exprs$dim2
    lines.sub <- concaveman::concaveman(cbind(xcoord, ycoord),
                                        concavity = concavity,
                                        length_threshold = length.threshold)
    
    concave.hulls <- rbind(concave.hulls, cbind(lines.sub, k = k))
  }
  
  concave.hulls <- data.frame(concave.hulls)
  colnames(concave.hulls) <- c("dim1", "dim2", "clusters")
  
  concave.hulls$dim1 <- as.numeric(concave.hulls$dim1)
  concave.hulls$dim2 <- as.numeric(concave.hulls$dim2)
  
  return(concave.hulls)
}

# @title Internal - Computes the number of cells for each cluster
#
# @description This function is used internally to computes the number of cells associated to each cluster.
#
# @param proj a data.frame providing the manifold representation with two columns
# @param clusters a character vector providing the cluster to analyse the associated with cell cluster
# @param samples a character vector providing for each cell the associated biological sample
#
# @return a data.frame containing the numbers of cells associated for each cluster for each sample (rownames = clusters / colnames = samples)
#
computeCellCounts <- function(proj,
                              clusters,
                              samples) {
  
  data <- cbind(proj, clusters, samples)
  
  cell.count <- plyr::ddply(data, c("clusters", "samples"), nrow)
  colnames(cell.count) <- c("clusters", "samples", "number")
  cell.count <- reshape::cast(cell.count, clusters ~ samples, value = "number")
  cell.count[is.na(cell.count)] <- 0
  rownames(cell.count) <- cell.count$clusters
  cell.count$clusters <- NULL
  cell.count <- data.frame(cell.count)
  
  return(cell.count)
}

# @title Internal - Computes the abundances for each cell cluster
#
# @description This function is used internally to computes the abundance of each cluster for each sample
#
# @param count a data.frame providing the numbers of cells associated to each cluster for each sample
#
# @return a data.frame containing the abundance of cells to each clusters for each sample
#
computeClusterAbundances <- function(count) {
  
  matrix.abundance <- apply(count, 2, function(df) {
    df / sum(df) }) * 100
  matrix.abundance <- data.frame(matrix.abundance)
  
  return(matrix.abundance)
}

#' @title Creates metaclusters
#'
#' @description This function aims to gathered multiple cell clusters to a large cell cluster
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to gather
#' @param metacluster.name a character value containing the name of the metacluster to create 
#'
#' @return a Celldata object
#'
#' @export
#'
createMetaclusters <- function(Celldata, 
                               clusters,
                               metacluster.name) {
  
  checkmate::qassert(clusters, "S+")
  checkmate::qassert(metacluster.name, "S1")
  
  identify.clusters <- Celldata@identify.clusters
  
  matrix.cell.count <- Celldata@matrix.cell.count
  matrix.cell.count <- cbind(matrix.cell.count, "clusters" = rownames(matrix.cell.count))
  
  for (i in clusters) {
    identify.clusters[identify.clusters == i] <- metacluster.name
    matrix.cell.count$clusters[matrix.cell.count$clusters == i] <- metacluster.name
  }
  
  matrix.cell.count <- plyr::ddply(matrix.cell.count, "clusters", function(x){
    x$clusters = NULL
    apply(x,2,sum)
  })
  
  rownames(matrix.cell.count) <- matrix.cell.count$clusters
  matrix.cell.count$clusters  <- NULL
  
  Celldata@identify.clusters <- identify.clusters
  Celldata@matrix.cell.count <- matrix.cell.count
  Celldata@matrix.abundance  <- computeClusterAbundances(count = matrix.cell.count)
  
  return(Celldata)
}

#' @title Deletes cluster
#'
#' @description This function aims to delete a set of cell clusters from this analysis
#'
#' @param Celldata a Celldata object
#' @param clusters a character vector containing the identifiers of the clusters to delete
#'
#' @return a Celldata object
#'
#' @export
#'
deleteClusters <- function(Celldata, 
                           clusters) {
  
  checkmate::qassert(clusters, "S+")
  
  identify.clusters <- Celldata@identify.clusters
  
  matrix.cell.count <- Celldata@matrix.cell.count
  matrix.cell.count <- cbind(matrix.cell.count, "clusters" = rownames(matrix.cell.count))
  
  for(i in clusters) {
    identify.clusters[identify.clusters == i] <- "todel"
    matrix.cell.count$clusters[matrix.cell.count$clusters == i] <- "todel"
  }
  
  matrix.cell.count <- plyr::ddply(matrix.cell.count, "clusters", function(x){
    x$clusters = NULL
    apply(x,2,sum)
  })
  
  matrix.cell.count$todel = NULL
  
  rownames(matrix.cell.count) <- matrix.cell.count$clusters
  matrix.cell.count$clusters  <- NULL
  
  Celldata@identify.clusters <- identify.clusters
  Celldata@matrix.cell.count <- matrix.cell.count
  Celldata@matrix.abundance  <- computeClusterAbundances(count = matrix.cell.count)
  
  return(Celldata)
}
