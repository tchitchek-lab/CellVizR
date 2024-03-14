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
    markernames[grepl("FS|SS", markernames)] <- names(markernames)[grepl("FS|SS", markernames)]
    
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
                           probs = c(0.05, 0.95)) {

  checkmate::qassert(files, "S*")
  checkmate::qassert(probs, "N2")

  table.check.min <- data.frame()
  table.check.max <- data.frame()
  table.check.median = data.frame()
  
  for (file in files) {

    fcs <- flowCore::read.FCS(file)
    trans.logicle <- flowCore::logicleTransform(w = 0.5, t = 262144, m = 4.5)
    marker.trans <- flowCore::colnames(fcs)
    trans <- flowCore::transformList(marker.trans, trans.logicle)
    fcs <- flowCore::transform(fcs, trans)
    data.expr <- flowCore::exprs(fcs)
    data.expr <- data.expr[, colnames(data.expr) %in%
                             names(flowCore::markernames(fcs))]

    ranges <- apply(data.expr, 2, stats::quantile, probs = probs)
    medians = apply(data.expr, 2, stats::median)
    
    table.check.min <- rbind(table.check.min, ranges[1, ])
    table.check.max <- rbind(table.check.max, ranges[2, ])
    table.check.median = rbind(table.check.median, medians)
  }

  markernames <- flowCore::markernames(fcs)

  colnames(table.check.min) <-  markernames
  colnames(table.check.max) <-  markernames
  colnames(table.check.median) <-  markernames
  
  rownames(table.check.min) <- gsub(".fcs", "", basename(files))
  rownames(table.check.max) <- gsub(".fcs", "", basename(files))
  rownames(table.check.median) <- gsub(".fcs", "", basename(files))
  
  res <- list(table.check.min, table.check.max, table.check.median)

  names(res)[1] <- paste0("quantiles_", probs[1])
  names(res)[2] <- paste0("quantiles_", probs[2])
  names(res)[3] <- "medians"
  
  
  
  graphData = NULL
  
  
  for(d in 1:length(res))
  {
    
    
    currentTable = data.frame(t(res[[d]]))
    
    for(s in 1:ncol(currentTable))
    {
      currentElement = data.frame(matrix(0, ncol = 4, nrow = nrow(currentTable)))
      colnames(currentElement) = c("value", "sample", "marker", "data")
      
      
      currentElement$value = currentTable[, s]
      currentElement$sample = colnames(currentTable)[s]
      currentElement$marker = rownames(currentTable)
      currentElement$data = names(res)[d]
      
      graphData = rbind(graphData, currentElement)
    }
    
    
    
  }
  
  # graphData = graphData[order(graphData$marker), ]
  
  plot <- ggplot2::ggplot(graphData, ggplot2::aes(x = value, y = marker, color = data, shape = data)) +
    
    ggplot2::geom_point(size = 1.5, alpha = 1) +
    ggplot2::geom_jitter(size = 1.5) +
    ggplot2::geom_boxplot(color = "grey20", linewidth = 0.75, alpha = 0,  position = ggplot2::position_dodge(width = 0)) +
    ggplot2::scale_color_discrete(labels = c("Median", "5th percentile", "95th percentile")) +
    ggplot2::scale_shape_discrete(labels = c("Median", "5th percentile", "95th percentile")) +
    
    ggplot2::labs(title = "Individualized marker ranges", x = "MFI/MSI", y = NULL) +
    ggplot2::theme_classic() +
    
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5)) +
  
    ggplot2::guides(shape = ggplot2::guide_legend(title = NULL), col = ggplot2::guide_legend(title = NULL))
  
  
  
  
  return(plot)
}