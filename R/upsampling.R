#' @title Performs the upsampling of downsampled events
#'
#' @description This function aims to perform the upsampling of downsampled events events based on an existing Celldata object and existing cell events stored in tab-separated or FCS files.
#'
#' Importantly, the identification of cell clusters must have been performed prior to this operation.
#'
#' @param Celldata a Celldata object
#' @param files a character vector providing the path of the tab-separated or FCS files
#' @param transform a character value containing the type of the transformation to apply. Possible values are: 'logicle', 'arcsinh', 'logarithmic' or 'none'
#'
#' @return a S4 object of class 'Celldata'
#' 
#' @export
#'
performUpsampling <- function(Celldata,
                              files,
                              transform = c("logicle", "arcsinh",
                                            "logarithmic", "none")) {

  checkmate::qassert(files, "S*")

  sample_files <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(files))
  files <- files[sample_files %in% unique(Celldata@samples)]
  Celldata.reload <- import(files, transform = transform)
  downsampled.exp <- Celldata@matrix.expression
  newmarkersnames <- colnames(downsampled.exp)
  colnames(downsampled.exp) <- Celldata@raw.markers

  reload.exp <- Celldata.reload@matrix.expression
  reload.exp <- reload.exp[, colnames(reload.exp) %in% colnames(downsampled.exp)]

  chk.downsampled <- apply(downsampled.exp, 1, sum)
  chk.reload <- apply(reload.exp, 1, sum)

  upsampled.exp <- reload.exp[!chk.reload %in% chk.downsampled, ]
  upsampled.samples <- Celldata.reload@samples[!chk.reload %in% chk.downsampled]
  
  downsampled.exp$cluster   <- Celldata@identify.clusters
  downsampled.centers <- plyr::ddply(downsampled.exp, "cluster", function(x) {
    x$cluster <- NULL
    centers <- apply(x, 2, stats::median, rm.na = TRUE)
    return(centers)
  })
  downsampled.centers$cluster <- NULL
  downsampled.exp$cluster <- NULL

  knn <- FNN::knnx.index(downsampled.centers, upsampled.exp,
                         k = 1, algorithm = "kd_tree")

  Celldata@matrix.expression <- rbind(downsampled.exp, upsampled.exp)
  colnames(Celldata@matrix.expression) <- newmarkersnames

  matrix.expression.r <- Celldata.reload@matrix.expression.r
  Celldata@matrix.expression.r <- matrix.expression.r[, colnames(matrix.expression.r) %in% colnames(downsampled.exp)]
  Celldata@samples <- c(Celldata@samples, upsampled.samples)
  Celldata@identify.clusters <- c(Celldata@identify.clusters, knn)

  message("computing cell cluster count matrix...")
  Celldata@matrix.cell.count <- computeCellCounts(proj = Celldata@matrix.expression,
                                                  clusters = Celldata@identify.clusters,
                                                  samples = Celldata@samples)

  message("computing cell cluster abundance matrix...")
  Celldata@matrix.abundance <- computeClusterAbundances(count = Celldata@matrix.cell.count)

  validObject(Celldata)
  return(Celldata)
  
}
