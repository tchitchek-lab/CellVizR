#' @title UMAPdata class definition
#'
#' @description The UMAPdata object is a S4 object containing all cytometry expressions.
#'
#' @slot samples a character vector containing the names of the biological samples
#' @slot raw.markers a character vector containing the names of the raw markers
#' @slot matrix.expression.r a data.frame containing the raw marker expressions of each cell
#' @slot matrix.expression a data.frame containing the marker expressions of each cell
#' @slot manifold a data.frame containing the manifold coordinates
#' @slot manifold.params a list containing the parameters used for manifold creation
#' @slot identify.clusters a vector containing the identifiers of cell clusters
#' @slot identify.clusters.params a vector containing the parameters used for the identification of the cell clusters
#' @slot concave.hulls a data.frame containing the coordinates of the cell cluster of the concave hulls for each cluster
#' @slot matrix.cell.count a data.frame containing the number of cells associated to each cluster for each sample
#' @slot matrix.abundance a data.frame containing the percentage of cells associated to each cluster for each sample
#' @slot statistic a data.frame containing the statistics of cell clusters
#' @slot metadata a data.frame containing the metadata associated to each sample
#'
#' @name UMAPdata-class
#' @rdname UMAPdata-class
#' @exportClass UMAPdata
#'

UMAPdata <- methods::setClass("UMAPdata",
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
                                # if (nrow(matrix.cell.count) != length(unique(samples))) {
                                #        stop("...")
                                # }
                                # if (nrow(matrix.abundance) != length(unique(samples))) {
                                #        stop("...")
                                # }
                                return(TRUE)
                              }
)