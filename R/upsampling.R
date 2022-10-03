#' @title Performs the upsampling of downsampled events
#'
#' @description This function aims to perform upsample downsampled events based on an existing UMAPdata object and existing cell events stored in tab-separated or FCS files
#'
#' Importantly, the identification of cell clusters must have been performed prior to this operation
#'
#' @param UMAPdata a UMAPdata object
#' @param files a character vector providing the path of the tab-separated or FCS files
#' @param transform a character value containing the type of the transformation to apply. Possible values are: 'logicle', 'arcsinh', 'logarithmic' or 'none'
#'
#' @return a S4 object of class 'UMAPdata'
performUpsampling <- function(UMAPdata,
                              files,
                              transform = c("logicle", "arcsinh",
                                            "logarithmic", "none")) {

  checkmate::qassert(files, "S*")

  sample_files <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(files))
  files <- files[sample_files %in% unique(UMAPdata@samples)]
  UMAPdata.reload <- import(files, transform = transform)
  downsampled.exp <- UMAPdata@matrix.expression
  newmarkersnames <- colnames(downsampled.exp)
  colnames(downsampled.exp) <- UMAPdata@raw.markers

  reload.exp <- UMAPdata.reload@matrix.expression
  reload.exp <- reload.exp[, colnames(reload.exp) %in% colnames(downsampled.exp)]

  chk.downsampled <- apply(downsampled.exp, 1, sum)
  chk.reload <- apply(reload.exp, 1, sum)

  upsampled.exp <- reload.exp[!chk.reload %in% chk.downsampled, ]
  upsampled.samples <- UMAPdata.reload@samples[!chk.reload %in% chk.downsampled]
  
  downsampled.exp$cluster   <- UMAPdata@identify.clusters
  downsampled.centers <- plyr::ddply(downsampled.exp, "cluster", function(x) {
    x$cluster <- NULL
    centers <- apply(x, 2, stats::median, rm.na = TRUE)
    return(centers)
  })
  downsampled.centers$cluster <- NULL
  downsampled.exp$cluster <- NULL

  knn <- FNN::knnx.index(downsampled.centers, upsampled.exp,
                         k = 1, algorithm = "kd_tree")

  UMAPdata@matrix.expression <- rbind(downsampled.exp, upsampled.exp)
  colnames(UMAPdata@matrix.expression) <- newmarkersnames

  matrix.expression.r <- UMAPdata.reload@matrix.expression.r
  UMAPdata@matrix.expression.r <- matrix.expression.r[, colnames(matrix.expression.r) %in% colnames(downsampled.exp)]
  UMAPdata@samples <- c(UMAPdata@samples, upsampled.samples)
  UMAPdata@identify.clusters <- c(UMAPdata@identify.clusters, knn)

  message("computing cell cluster count matrix...")
  UMAPdata@matrix.cell.count <- computeCellCounts(proj = UMAPdata@matrix.expression,
                                                  clusters = UMAPdata@identify.clusters,
                                                  samples = UMAPdata@samples)

  message("computing cell cluster abundance matrix...")
  UMAPdata@matrix.abundance <- computeClusterAbundances(count = UMAPdata@matrix.cell.count)

  validObject(UMAPdata)
  return(UMAPdata)
}
