#' @title Identify cell cluster of having similar marker expressions
#'
#' @description This function aims to identify cell clusters, which are groups of cells having similar expressions for selected markers, using different unsupervised clustering methods.
#' 
#' Several clustering method are available such as kmeans, kmedian, clara, DBSCAN, HDBSCAN and SOM.
#' The cell clustering can be performed on the manifold representation or based on marker expression. 
#' 
#' @details 
#' For each identify cell cluster, the boundaries of cells belonging to this cluster are delineated using a concave hull
#'
#' @param UMAPdata a UMAPdata object
#' @param space a character value containing the space of clustering method to use. Possible values are: 'manifold' or 'markers'
#' @param method a character value containing the name of the clustering method to use. Possible values are: 'kmeans', 'kmedian', 'clara', 'DBSCAN' and 'SOM'
# @param k a numeric value providing the number of clusters to identify
# @param nstart a numeric value providing the number of random sets that should be chosen (please refer to the function 'kmeans' of the 'stats' package)
#' @param concavity a numeric value providing a relative measure of concavity for the computation of the concave hulls (please refer to the function 'concaveman' of the 'concaveman' package)
#' @param length.threshold a numeric value providing a threshold of the segment length for the computation of the concave hulls (please refer to the function 'concaveman' of the concaveman package)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param ... Other arguments passed on to methods
#'
#' @return a S4 object of class 'UMAPdata'
#'
#' @export
#'
identifyClusters <- function(UMAPdata,
                             space = c("manifold", "markers"),
                             method = c("kmeans", "kmedian", "clara", "DBSCAN", "SOM"),
                             concavity = 2,
                             length.threshold = 0,
                             seed = 42,
                             ...) {
  
  space <- match.arg(space)
  method <- match.arg(method)
  
  checkmate::qassert(space, "S1")
  checkmate::qassert(method, "S1")
  checkmate::qassert(concavity, "N1")
  checkmate::qassert(length.threshold, "N1")
  checkmate::qassert(seed, "N1")
  
  proj  <- UMAPdata@manifold
  exprs <- UMAPdata@matrix.expression
  
  message(paste0("Identifying cell clusters..."))
  
  if(space=="manifold"){
    data = proj
  }else if (space == "markers") {
    data = exprs
  }
  
  set.seed(seed)
  switch(method,
         kmeans = {
           clusters <- stats::kmeans(data, ...)$cluster
         },
         kmedian = {
           clusters = Gmedian::kGmedian(data, ...)$cluster[,1]
         },
         clara = {
           clusters <- cluster::clara(data, ...)$clustering
         },
         DBSCAN = {
           clusters <- dbscan::dbscan(data, ...)$cluster 
         },
         SOM = {
           clusters <- kohonen::som(as.matrix(data), ...)$unit.classif # grid(5*5, hexagonal)
         })
  
  clusters <- as.vector(clusters)
  clusters <- as.character(clusters)
  
  UMAPdata@identify.clusters <- clusters
  
  identify.Clusters.params <- list(...,
                                   concavity = concavity,
                                   length.threshold = length.threshold)
  UMAPdata@identify.clusters.params <- identify.Clusters.params
  
  if (space == "manifold") {
    message(paste0("computing cell clusters boundaries..."))
    UMAPdata@concave.hulls <- computeConcaveHulls(proj = proj,
                                                  clusters = clusters,
                                                  concavity = concavity,
                                                  length.threshold = length.threshold)
  }
  
  samples <- UMAPdata@samples
  
  message(paste0("computing cell cluster count matrix..."))
  UMAPdata@matrix.cell.count <- computeCellCounts(proj = data,
                                                  samples = samples,
                                                  clusters = clusters)
  
  message(paste0("computing cell cluster abundance matrix..."))
  count <- UMAPdata@matrix.cell.count
  UMAPdata@matrix.abundance <- computeClusterAbundances(count = count)
  
  
  return(UMAPdata)
  
}

#' @title Internal - Computes the concave hulls for identified cell clusters
#'
#' @description This function is used internally to computes the concave hulls each cell cluster.
#' Each concave hulls correspond the boundaries of the cells cluster with manifold representation.
#'
#' @param proj a data.frame providing the manifold representation
#' @param clusters a character vector providing id of cell clusters for which the concave hulls must be computed
#' @param concavity a numeric value providing a relative measure of concavity (please refer to the function 'concaveman' of the 'concaveman' package)
#' @param length.threshold a numeric value providing a threshold for the segment length (please refer to the function 'concaveman' of the concaveman package)
#'
#' @return a data.frame containing the concave hulls for each cluster (dim1, dim2, clusters)
#'
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

#' @title Internal - Computes the number of cells for each cluster
#'
#' @description This function is used internally to computes the number of cells for each cluster.
#'
#' @param proj a data.frame providing the manifold representation with two columns
#' @param clusters a character vector providing the cluster to analyse the associated with cell cluster
#' @param samples a character vector providing for each cell the associated biological sample
#'
#' @return a data.frame containing the numbers of cells associated for each cluster for each sample (rownames = clusters / colnames = samples)
#'
computeCellCounts <- function(proj,
                              clusters,
                              samples) {
  
  
  data <- cbind(proj, clusters, samples)
  
  cell.count <- plyr::ddply(data, c("clusters", "samples"), nrow)
  colnames(cell.count) <- c("clusters", "samples", "number")
  cell.count <- reshape::cast(cell.count, clusters~samples, value = "number")
  cell.count[is.na(cell.count)] <- 0
  rownames(cell.count) <- cell.count$clusters
  cell.count$clusters <- NULL
  cell.count <- data.frame(cell.count)
  
  return(cell.count)
}

#' @title Internal - Computes the abundances for each cell cluster
#'
#' @description This function is used internally to computes the abundance of each cluster for each sample.
#'
#' @param count a data.frame providing the numbers of cells associated to each cluster for each sample
#'
#' @return a data.frame containing the abundance of cells to each clusters for each sample
#'
computeClusterAbundances <- function(count) {
  
  matrix.abundance <- apply(count, 2, function(df){df / sum(df)}) * 100
  matrix.abundance <- data.frame(matrix.abundance)
  
  return(matrix.abundance)
  
}
