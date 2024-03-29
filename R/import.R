#' @title Imports cell expression profiles from TSV or FCS files
#'
#' @description This function aims to import acquired cell events from cytometric profiling into a Celldata object.
#'
#' Input files can be tab-separated or FCS files.
#' Different transformations can be applied such as logicle, arcsinh or logarithmic.
#' Importantly, a downsampling of cell events can be performed using uniformly-based or density-based random selections.
#' Cell marker having technical or biological biaises can be excluded during the import.
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#' @param filetype a character vector specifying the format of the loaded files. By default, FCS is used
#' @param transform a character value containing the type of the transformation to apply. Possible values are: 'logicle', 'arcsinh', 'logarithmic' or 'none'
#' @param d.method a character value containing the type of the downsampling to apply. Possible values are: 'none' or 'uniform'
#' @param parameters.method a list value containing the parameters to use for downsampling
#' @param exclude.markers a character vector providing the marker names to be excluded during the import
#' @param seed a numeric value providing the random seed to use during stochastic operations
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#' @import methods
import <- function(files,
                   filetype = "fcs",
                   transform = c("logicle", "arcsinh", "logarithmic", "none"),
                   d.method = c("none", "uniform"),
                   parameters.method = list("target.number" = NULL, "target.percent" = 0.1),
                   exclude.markers = NULL,
                   seed = 42) {

  transform <- match.arg(transform)
  d.method <- match.arg(d.method)

  checkmate::qassert(filetype, "S1")
  checkmate::qassert(transform, "S1")
  if (d.method == "uniform") {
    checkmate::qassert(parameters.method, "L+")
  }
  checkmate::qassert(d.method, "S1")
  checkmate::qassert(exclude.markers, c("0", "S*"))
  checkmate::qassert(seed, "N1")

  exprs.r <- data.frame()
  exprs <- data.frame()
  samples <- c()

  message("Transformation method is: ", transform)
  if (d.method == "uniform") {
    message("Downsampling method is: ", d.method, " and parameters method are: ",
            paste0(names(parameters.method), " = ", parameters.method, collapse = ", "))
    message()
  } else {
    message("Downsampling method is: ", d.method)
    message()
  }

  for (file in files) {
    message("importing ", basename(file), " file")

    if (filetype == "fcs") {
      fcs <- flowCore::read.FCS(file)
    } else {
      exprs <- utils::read.delim(file, sep = "\t")
      exprs <- exprs[, !colnames(exprs) %in% exclude.markers]
      fcs <- suppressWarnings(createFlowframe(exprs))
    }

    exprs.raw <- flowCore::exprs(fcs)
    exprs.raw <- exprs.raw[, colnames(exprs.raw) %in% names(flowCore::markernames(fcs))]

    switch(transform,
           logicle = {
             trans.logicle <- flowCore::logicleTransform(w = 0.5, t = 262144, m = 4.5)
             marker.trans <- flowCore::colnames(fcs)
             trans <- flowCore::transformList(marker.trans, trans.logicle)
             fcs <- flowCore::transform(fcs, trans)
           },
           arcsinh = {
             trans.arcsinh <- flowCore::arcsinhTransform(a = 0, b = 0.2)
             marker.trans <- flowCore::colnames(fcs)
             trans <- flowCore::transformList(marker.trans, trans.arcsinh)
             fcs <- flowCore::transform(fcs, trans)
           },
           logarithmic = {
             trans.log <- flowCore::logTransform(logbase = 10, r = 1, d = 1)
             marker.trans <- flowCore::colnames(fcs)
             trans <- flowCore::transformList(marker.trans, trans.log)
             fcs <- flowCore::transform(fcs, trans)
           },
           none = {
             fcs <- fcs
           })

    exprs.sub <- flowCore::exprs(fcs)
    exprs.sub <- exprs.sub[, colnames(exprs.sub) %in% names(flowCore::markernames(fcs))]

    if (d.method == "none") {
      exprs.sub <- exprs.sub

    } else if (d.method == "uniform") {
      sampled <- downsamplingUniform(file = file,
                                     exprs.raw = exprs.raw,
                                     target.percent = parameters.method$target.percent,
                                     target.number = parameters.method$target.number,
                                     seed = seed)

      exprs.raw <- exprs.raw[sampled, ]
      exprs.sub <- exprs.sub[sampled, ]

    }

    if (file == files[1]) {
      markernames <- flowCore::markernames(fcs)
      markernames[grepl("FS|SS", markernames)] <- names(markernames)[grepl("FS|SS", markernames)]

      colnames(exprs.raw) <- markernames
      colnames(exprs.sub) <- markernames

    } else {
      colnames(exprs.raw) <- colnames(exprs.r)
      colnames(exprs.sub) <- colnames(exprs)
    }

    exprs.r <- rbind(exprs.r, exprs.raw)
    exprs <- rbind(exprs, exprs.sub)

    sample <- gsub(".fcs", "", basename(file))
    samples <- c(samples, rep(sample, nrow(exprs.sub)))

  }

  if (!is.null(exclude.markers)) {
    exprs.r <- exprs.r[, !colnames(exprs.r) %in% exclude.markers]
    exprs   <- exprs[, !colnames(exprs) %in% exclude.markers]
  }

  res <- methods::new("Celldata",
                      matrix.expression.r = exprs.r,
                      matrix.expression = exprs,
                      samples = samples,
                      raw.markers = colnames(exprs),
                      matrix.abundance = data.frame())

  return(res)
}

#' @title Imports cell expression profiles from MTX files
#'
#' @description This function aims to import acquired cell events from single-cell transcriptomic profiling into a Celldata object.
#'
#' @param count a character vector specifying the path of the count mtx file containing the count values
#' @param cells a character vector specifying the path of the mtx file containing the cell ids
#' @param features a character vector specifying the path of the mtx file containing the gene ids
#' @param meta a character specifying the name of the tsv metadata file to load
#' @param sample.col a character specifying the name of the colmun to use as biological sample in the metadata file
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#' @import methods 
importMTX <- function(count,
                      cells,
                      features,
					  meta,
					  sample.col) {
		
	expression_matrix <- Seurat::ReadMtx(mtx      = count,
							             cells    = cells,
							             features = features)
							 
	seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix)
	counts        <- Seurat::GetAssayData(object = seurat_object[["RNA"]], slot = "counts")
	counts        <- t(as.matrix(x = log(counts + 1)))
	counts        <- counts[order(rownames(counts)),]
	counts        <- data.frame(counts)
	
	meta                <- utils::read.delim(meta,header=TRUE)
	meta                <- meta[meta$id %in% rownames(counts),]
	samples             <- meta[,sample.col,drop=TRUE]
	
	res <- methods::new("Celldata",
                      matrix.expression.r = counts,
                      matrix.expression = counts,
                      samples = samples,
                      raw.markers = colnames(counts),
                      matrix.abundance = data.frame())

	return(res)
	
}

# @title Internal - Computes the downsampling uniformly-based
#
# @description This function aims to perform a downsampling uniformly-based cells
#
# @param file a character vector specifying the path of the tab-separated or FCS files to load
# @param exprs.raw a data.frame providing the raw marker expressions
# @param target.percent a numeric value providing the percentage of cells to downsample for each sample
# @param target.number a numeric value providing the number of cells to downsample for each sample
# @param seed a numeric value providing the random seed to use during stochastic operations
#
# @return an index vector for the downsampling
#
# @export
#
downsamplingUniform <- function(file,
                                exprs.raw,
                                target.percent,
                                target.number,
                                seed) {

  do.call("set.seed", list(seed))

  if (!is.null(target.percent) && is.null(target.number)) {
    target.number <- round(target.percent * nrow(exprs.raw))
    message("Number of cells: ", target.number, " for ", basename(file))
    sampled <- sample(nrow(exprs.raw), min(target.number, nrow(exprs.raw)))
  } else if (is.null(target.percent) && !is.null(target.number)) {
    sampled <- sample(nrow(exprs.raw), min(target.number, nrow(exprs.raw)))
  } else if (is.null(target.percent) && is.null(target.number)) {
    stop("Parameters are null")
  } else if (!is.null(target.percent) && !is.null(target.number)) {
    stop("Specify only one parameter")
  }
  return(sampled)
}

#' @title Renames markers within a Celldata object
#'
#' @description This function aims to rename cell markers stored within a Celldata object.
#'
#' This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
#'
#' @param Celldata a Celldata object
#' @param marker.names a character vector providing the new marker names to use
#'
#' @return a S4 object of class 'Celldata'
#'
#'@export
#'
renameMarkers <- function(Celldata,
                         marker.names) {

  if (length(marker.names) != length(unique(marker.names))) {
    stop("The same marker is used twice")
  }

  if (ncol(Celldata@matrix.expression) != length(marker.names)) {
    stop("The number of markers given is different from the initial number of markers")
  }

  colnames(Celldata@matrix.expression) <- marker.names
  validObject(Celldata)

  return(Celldata)
}

#' @title Assigns meta-information to biological samples
#'
#' @description This function aims to attach meta-information to biological samples.
#'
#' For each biological sample, the biological individual, the biological condition and the time point can be specified for subsequent analyses.
#'
#' @param Celldata a Celldata object
#' @param metadata a data.frame containing contextual information about the biological samples. This data.frame must have 3 columns specifying for each sample the associated individual (column named 'individual'), the biological condition (column named 'condition') and the time point (column named 'timepoint'). Rownames must correspond to biological samples imported within the Celldata object.
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#'
assignMetadata <- function(Celldata,
                          metadata) {

  checkmate::qassert(metadata, c("D2", "D3"))

  if (length(unique(colnames(metadata))) != length(colnames(metadata))) {
    stop("colnames are not unique")
  }

  if (!all(colnames(metadata) %in% c("individual", "condition", "timepoint"))) {
    stop("colnames must be 'individual', 'condition', or 'timepoint'")
  }

  if (!all(c("individual") %in% colnames(metadata))) {
    stop("colnames 'individual' is missing")
  }

  if (all(unique(Celldata@samples) != rownames(metadata))) {
    stop("sample names (rownames) of the metadata are inconcistent with the sample names stored in the Celldata object")
  }

  Celldata@metadata <- metadata

  validObject(Celldata)
  return(Celldata)
}
