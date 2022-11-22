# @title Internal - Computes the cell density within the manifold
#
# @description This function is used internally to compute the cell density.
#
# @param proj a data.frame providing the manifold representation with two columns
# @param n a numerical providing the number a grid points in each direction (please refer to the function 'kde2d' of the 'MASS' package)
#
# @return a numeric vector containing the computed cell density
#
computeCellDensities <- function(proj,
                                 n) {
  
  density <- MASS::kde2d(proj[, 1], proj[, 2], n)
  ix <- findInterval(proj[, 1], density$x)
  iy <- findInterval(proj[, 2], density$y)
  ii <- cbind(ix, iy)
  
  return(density$z[ii])
}

# @title Internal - Rescales numeric matrix for projection
#
# @description This functions is used internally to rescale expression values.
#
# @param proj a data.frame providing the manifold representation
# @param abs a boolean value specifying if data must be rescaled
# @param quant.low a numeric value providing the number first quantile
# @param quant.high a numeric value providing the number last quantile
#
# @return a numeric vector containing scale limits
#
abs.proj <- function(proj,
                     abs,
                     quant.low,
                     quant.high) {
  if (abs == FALSE) {
    quantiles <- stats::quantile(proj$value, probs = c(quant.low, quant.high))
    proj$value[proj$value < quantiles[1]] <- quantiles[1]
    proj$value[proj$value > quantiles[2]] <- quantiles[2]
    limits <- c(quantiles[1], quantiles[2])
  } else {
    proj$value[proj$value < quant.low] <- quant.low
    proj$value[proj$value > quant.high] <- quant.high
    limits <- c(quant.low, quant.high)
  }
  
  return(limits)
}

#' @title Select samples based on metadata information
#'
#' @description This function aims to select biological samples of interest based on provided metadata.
#'
#' @param Celldata a Celldata object
#' @param individual a character vector indicating the individual to select
#' @param condition a character vector indicating the biological condition to select
#' @param timepoint a character vector indicating the timepoint to select
#'
#' @return a character vector
#'
#' @export
#'
selectSamples <- function(Celldata,
                          individual = NULL,
                          condition = NULL,
                          timepoint = NULL) {
  
  metadata <- Celldata@metadata
  
  if(!is.null(individual)) {
    metadata <- metadata[metadata$individual %in% individual, ]
  }
  if(!is.null(condition)) {
    metadata <- metadata[metadata$condition %in% condition, ]
  }
  if(!is.null(timepoint)) {
    metadata <- metadata[metadata$timepoint %in% timepoint, ]
  }
  
  samples <- rownames(metadata)
  
  return(samples)
} 
