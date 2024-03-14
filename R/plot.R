#' @title Plots graphics for all Celldata objects
#'
#' @description Generates a graphical representation for a Celldata object.
#' The displayed representation depends on the current analysis status of the Celldata object.
#'
#' - If the manifold has not been calculated, then the number of cells per sample will be displayed
#'
#' - If the manifold has been calculated but not the clustering, then the manifold representation will be displayed
#'
#' - If the manifold and clustering have been calculated, then a heatmap of marker expressions will be displayed
#'
#' @param x a Celldata object
#'
#' @return a ggplot2 object
#'
#' @name plot
#' @rdname plot-methods

setMethod("plot", "Celldata",
          function(x) {
            if (length(x@manifold) == 0 && length(x@identify.clusters) == 0) {
              plot(plotCellCounts(x))
            } else if (length(x@manifold) != 0 && length(x@identify.clusters) == 0) {
              plot(plotManifold(x))
            } else if (length(x@manifold) != 0 && length(x@identify.clusters) != 0) {
              plot(plotHmExpressions(x, markers = colnames(x@matrix.expression))$hm)
            }
          })
