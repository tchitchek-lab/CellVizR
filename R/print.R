#' @title Prints information for a given UMAPdata object
#'
#' @description Prints a preview for a UMAPdata object.
#'  
#' @param x a UMAPdata object 
#'  
#' @return none
#'  
#' @name print
#' @rdname print-methods 
NULL
#' @rdname print-methods

setMethod("print", "UMAPdata", 
          function(x){
            cat("Object class: UMAPdata\n")
            cat(paste0("Numbers of markers: ",length(x@matrix.expression[,-1])))
            cat("\nMarkers: ")
            cat(paste0(colnames(x@matrix.expression), collapse = ", "))
            cat("\n")
            cat(paste0("Numbers of samples: ",length(unique(x@samples))))
            cat("\nSamples: ")
            cat(paste0(unique(x@samples)[seq(1,3)], collapse = ", "), "...") 
            cat("\n")
            cat("Numbers of cells: ")
            cat(formatC(nrow(x@matrix.expression), big.mark = ","))
            cat("\n")
            cat("- Metadata: ")
            cat("\n")
            if(all(c("condition") %in% colnames(x@metadata))) { 
              cat(paste0("Conditions: ", paste0(names(table(x@metadata$condition )), "=", table(x@metadata$condition) , collapse = ", ")))
              cat("\n")
            }
            if(all(c("timepoint") %in% colnames(x@metadata))) { 
              cat(paste0("Timepoint: ", paste0(names(table(x@metadata$timepoint)), "=", table(x@metadata$timepoint) , collapse = ", ")))
              cat("\n") 
            }
            cat("- Manifold")
            cat("\n")
            if (length(x@manifold) == 0) {
              cat("No manifold")
            } else {
              # cat(paste0("Manifold markers are: ", paste0(unlist(manifold.markers), collapse = ", ")))
              cat(paste0("Parameters: ", paste0(names(x@manifold.params), "=", x@manifold.params, collapse = ", " )))
            }
            cat("\n")
            cat("- Clustering")
            cat("\n")
            if (length(x@identify.clusters) == 0) {
              cat("No clustering")
            } else {
              cat(paste0("Numbers of clusters: ", length(unique(x@identify.clusters))))
              cat("\n")
              cat(paste0("Parameters: ", paste0(names(x@identify.clusters.params), "=", x@identify.clusters.params, collapse = ", " )))
              cat("\n")
            }
            
          })

#' @title Prints information for a UMAPdata objects 
#' 
#' @description Shows a preview for a UMAPdata object.
#' 
#' @param object a UMAPdata object 
#' 
#' @return none
#' 
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods

setMethod("show", "UMAPdata", 
          definition = function(object){print(object)}
)