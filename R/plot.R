#' @title Plots graphics for all UMAPdata objects 
#'
#' @description Generates a graphical representation for a UMAPdata object.
#' The displayed representation depends on the current analysis status of the UMAPdata object. 
#'
#' - If the manifold has not been calculated, then the number of cells per sample will be displayed.
#' 
#' - If the manifold has been calculated but not the clustering, then the manifold representation will be displayed. 
#' 
#' - If the manifold and clustering have been calculated, then a heatmap of marker expressions will be displayed.  
#'  
#' @param x a UMAPdata object 
#'  
#' @return a ggplot2 object
#'  
#' @name plot
#' @rdname plot-methods 
NULL
#' @rdname plot-methods

setMethod("plot", "UMAPdata", 
          function(x){
            if (length(x@manifold) == 0 & length(x@identify.clusters) == 0) {
              plot(plotCellCounts(x))
            } else if (length(x@manifold) != 0 & length(x@identify.clusters) == 0) {
              plot(plotManifold(x))
            } else if (length(x@manifold) != 0 & length(x@identify.clusters) != 0) {
              plot(plotHmExpressions(x, markers = colnames(x@matrix.expression))$hm)
            }
          })


