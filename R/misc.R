# @title Internal - Computes the cell density within the manifold
#
# @description This function is used internally to compute the cell density
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
# @description This functions is used internally xxx
#
# @param proj a data.frame providing the manifold representation
# @param abs a boolean value providing xx
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
