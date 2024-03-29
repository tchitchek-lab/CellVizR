#' @title Prints information for a given Celldata object
#'
#' @description Prints a preview for a Celldata object.
#'
#' @param x a Celldata object
#'
#' @return none
#'
#' @name print
#' @rdname print-methods
 
setMethod("print", "Celldata",
          function(x) {
            cat("Object class: Celldata\n")
            cat(paste0("Numbers of samples: ", length(unique(x@samples))))
            cat("\n- Samples: ")
            if(length(unique(x@samples))>10){
				cat(paste0(unique(x@samples)[seq(1, 10)], collapse = ", "), "...")
            }else{
				cat(paste0(unique(x@samples), collapse = ", "))
            }
			cat("\n")
            cat(paste0("Numbers of markers: ", formatC(length(x@matrix.expression[, -1]), big.mark = ","),"\n"))
            if(ncol(x@matrix.expression)<40){
				cat("- Markers: ")
				cat(paste0(colnames(x@matrix.expression), collapse = ", "))
				cat("\n")
			}
            cat("Numbers of cells: ")
            cat(formatC(nrow(x@matrix.expression), big.mark = ","))
            cat("\n")
            cat("- Metadata: ")
            cat("\n")
            if (all(c("condition") %in% colnames(x@metadata))) {
              cat(paste0("Conditions: ", paste0(names(table(x@metadata$condition)),
                                                "=", table(x@metadata$condition),
                                                collapse = ", ")))
              cat("\n")
            }
            if (all(c("timepoint") %in% colnames(x@metadata))) {
              cat(paste0("Timepoint: ", paste0(names(table(x@metadata$timepoint)),
                                               "=", table(x@metadata$timepoint),
                                               collapse = ", ")))
              cat("\n")
            }
            cat("- Manifold")
            cat("\n")
            if (length(x@manifold) == 0) {
              cat("No manifold")
              cat("\n")
            } else {
              cat(paste0("Parameters: ", paste0(names(x@manifold.params),
                                                "=", x@manifold.params,
                                                collapse = ", "),"\n"))
            }
            #cat("\n")
            cat("- Clustering")
            cat("\n")
            if (length(x@identify.clusters) == 0) {
              cat("No clustering")
              cat("\n")
            } else {
              cat(paste0("Numbers of clusters: ", length(unique(x@identify.clusters))))
              cat("\n")
              cat(paste0("Parameters: ", paste0(names(x@identify.clusters.params),
                                                "=", x@identify.clusters.params,
                                                collapse = ", ")))
              cat("\n")
            }
          })

#' @title Prints information for a Celldata objects
#'
#' @description Shows a preview for a Celldata object
#'
#' @param object a Celldata object
#'
#' @return none
#'
#' @name show
#' @rdname show-methods

setMethod("show", "Celldata",
          definition = function(object) {
            print(object)
          }
)
