#' @title Celldata class definition
#'
#' @description The Celldata object is a S4 object containing all single-cell information.
#'
#' @slot samples a character vector containing the names of the biological samples
#' @slot raw.markers a character vector containing the names of the raw markers
#' @slot matrix.expression.r a data.frame containing the raw marker expressions of each cell
#' @slot matrix.expression a data.frame containing the marker expressions of each cell
#' @slot manifold a data.frame containing the manifold coordinates
#' @slot manifold.params a list containing the parameters used for manifold generation
#' @slot identify.clusters a vector containing the names of the identified cell clusters
#' @slot identify.clusters.params a vector containing the parameters used for the identification of the cell clusters
#' @slot concave.hulls a data.frame containing the coordinates of the concave hulls of each cluster
#' @slot matrix.cell.count a data.frame containing the number of cells associated to each cluster for each sample
#' @slot matrix.abundance a data.frame containing the percentage of cells associated to each cluster for each sample
#' @slot statistic a data.frame containing the statistics of cell clusters
#' @slot metadata a data.frame containing the metadata associated to each sample
#'
#' @name Celldata-class
#' @rdname Celldata-class
#' @exportClass Celldata
#'
Celldata <- methods::setClass("Celldata",
                              slots = c(samples = "vector",
                                        raw.markers = "vector",
                                        matrix.expression.r = "data.frame",
                                        matrix.expression = "data.frame",
                                        manifold = "data.frame",
                                        manifold.params = "list",
                                        identify.clusters = "vector",
                                        identify.clusters.params = "list",
                                        concave.hulls = "data.frame",
                                        matrix.cell.count = "data.frame",
                                        matrix.abundance = "data.frame",
                                        statistic = "data.frame",
                                        metadata = "data.frame"),

                              validity = function(object) {
                                if (ncol(object@matrix.abundance) != 0 && ncol(object@matrix.abundance) != length(unique(object@samples))) {
                                  return("Error in Celldata object: The number of columns in the abundance matrix: ", ncol(object@matrix.abundance),
                                         ", is inconsistent with the number of samples: ", length(unique(object@samples)))
                                }
                                if (ncol(object@matrix.cell.count) != 0 && ncol(object@matrix.cell.count) != length(unique(object@samples))) {
                                  return("Error in Celldata object: The number of columns in the cell count matrix: ", ncol(object@matrix.cell.count),
                                         ", is inconsistent with the number of samples: ", length(unique(object@samples)))
                                }
                                if (nrow(object@matrix.abundance) != 0 && nrow(object@matrix.abundance) != length(unique(object@identify.clusters))) {
                                  return("Error in Celldata object: The number of row in the abundance matrix: ", nrow(object@matrix.abundance),
                                         ", is inconsistent with the number of clusters: ", length(unique(object@identify.clusters)))
                                }
                                if (nrow(object@matrix.cell.count) != 0 && nrow(object@matrix.cell.count) != length(unique(object@identify.clusters))) {
                                  return("Error in Celldata object: The number of row in the cell count matrix: ", nrow(object@matrix.cell.count),
                                         ", is inconsistent with the number of clusters: ", length(unique(object@identify.clusters)))
                                }
                                if (nrow(object@metadata) != 0 && nrow(object@metadata) != length(unique(object@samples))) {
                                  return("Error in Celldata object: The number of row in the metadata: ", nrow(object@metadata),
                                         ", is inconsistent with the number of samples: ", length(unique(object@samples)))
                                }
                                if (ncol(object@matrix.expression.r) != length(object@raw.markers)) {
                                  return("Error in Celldata object: The number of columns in the expression raw matrix: ", ncol(object@matrix.expression.r),
                                         ", is inconsistent with the number of markers: ", length(object@raw.markers))
                                }
                                if (ncol(object@matrix.expression) != length(object@raw.markers)) {
                                  return("Error in Celldata object: The number of columns in the expression matrix: ", ncol(object@matrix.expression),
                                         ", is inconsistent with the number of markers: ", length(object@raw.markers))
                                }
                                return(TRUE)
                              }
)
