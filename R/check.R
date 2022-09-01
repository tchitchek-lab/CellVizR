#' @title Verifies the consistency of the marker names within cell event files
#'
#' @description This function aims to check the consistency of marker names across multiple tab-separated or FCS files.
#'
#' Additionally, the number of cells associated to each sample is displayed.
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to check
#'
#' @return a data.frame containing the marker names and the associated number of cells for each sample (rownames = samples and colnames = markers)
#'
#' @export
#'
QCMarkerNames <- function(files) {
  
  checkmate::qassert(files, "S*")
  
  table.check <- data.frame()
  
  for (file in files) {
    fcs <- flowCore::read.FCS(file)
    data.expr <- flowCore::exprs(fcs)
    
    nbr.cells <- nrow(data.expr)
    markernames <- flowCore::markernames(fcs)
    
    table.check <- rbind(table.check, c(nbr.cells, markernames))
    
  }
  
  colnames(table.check) <- c("nb_cells", names(markernames))
  rownames(table.check) <- gsub(".fcs", "", basename(files))
  
  return(table.check)
}

#' @title Verifies the consistency of marker expressions integrity within cell event files
#'
#' @description This function aims to check the consistency of marker expressions ranges across multiple tab-separated or FCS files.
#'
#' The marker expressions ranges are calculated based on the user-defined quantiles.
#'
#' @param files a character vector specifying the path of the FCS files to verified
#' @param probs a numerical vector providing the quantiles used to define marker expressions ranges
#'
#' @return a list containing two data.frame for the lower and upper marker expression ranges (rownames = samples and colnames = markers)
#'
#' @export
#'
QCMarkerRanges <- function(files,
                           probs=c(0.05, 0.95)) {
  
  checkmate::qassert(files, "S*")
  checkmate::qassert(probs, "N2")
  
  table.check.min <- data.frame()
  table.check.max <- data.frame()
  
  for (file in files) {
    
    fcs <- flowCore::read.FCS(file)
    
    trans <- flowCore::estimateLogicle(fcs,
                                       channels = flowCore::colnames(fcs),
                                       m = 5.5)
    fcs <- flowCore::transform(fcs, trans)
    data.expr <- flowCore::exprs(fcs)
    data.expr <- data.expr[, colnames(data.expr) %in%
                             names(flowCore::markernames(fcs))]
    
    ranges <- apply(data.expr, 2, stats::quantile, probs = probs)
    
    table.check.min <- rbind(table.check.min, ranges[1, ])
    table.check.max <- rbind(table.check.max, ranges[2, ])
    
  }
  
  markernames <- flowCore::markernames(fcs)
  
  colnames(table.check.min) <-  markernames
  colnames(table.check.max) <-  markernames
  
  rownames(table.check.min) <- gsub(".fcs", "", basename(files))
  rownames(table.check.max) <- gsub(".fcs", "", basename(files))
  
  res <- list(table.check.min, table.check.max)
  
  names(res)[1] <- paste0("quantiles_", probs[1])
  names(res)[2] <- paste0("quantiles_", probs[2])
  
  return(res)
  
}
