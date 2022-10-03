#' @title Imports of cell expression profiles from TSV or FCS files
#'
#' @description This function aims to import acquired cell events into a UMAPdata object.
#'
#' Input files can be tab-separated or FCS files.
#' Different transformations can be applied such as logicle, arcsinh or logarithmic.
#' Importantly, a downsampling of cell events can be performed using uniformly-based or density-based random selections.
#' Cell marker having technical or biological biaises can be excluded during the import.
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#' @param filetype a character vector specifying the format of the loaded files. By default, FCS is used
#' @param transform a character value containing the type of the transformation to apply. Possible values are: 'logicle', 'arcsinh', 'logarithmic' or 'none'
#' @param d.method a character value containing the type of the downsampling to apply. Possible values are: 'none', 'uniform' or 'density'
#' @param parameters.method a list value containing the parameters to use for downsampling
#' @param exclude.markers a character vector providing the marker names to be excluded during the import
#' @param seed a numeric value providing the random seed to use during stochastic operations
#'
#' @return a S4 object of class 'UMAPdata'
#'
#' @export
#' @import methods
import <- function(files,
                   filetype = "fcs",
                   transform = c("logicle", "arcsinh", "logarithmic", "none"),
                   d.method = c("none", "uniform", "density"),
                   parameters.method = list("exclude.pctile" = 0.01, "target.pctile" = 0.05,
                                            "target.number" = NULL, "target.percent" = 0.1),
                   exclude.markers = NULL,
                   seed = 42) {

  transform <- match.arg(transform)
  d.method <- match.arg(d.method)

  checkmate::qassert(filetype, "S1")
  checkmate::qassert(transform, "S1")
  if (d.method == "uniform" || d.method == "density") {
    checkmate::qassert(parameters.method, "L+")
  }
  checkmate::qassert(d.method, "S1")
  checkmate::qassert(exclude.markers, c("0", "S*"))
  checkmate::qassert(seed, "N1")

  exprs.r <- data.frame()
  exprs <- data.frame()
  samples.r <- c()
  samples <- c()

  message("Transformation method is: ", transform)
  if (d.method == "uniform" || d.method == "density") {
    message("Downsampling method is: ", d.method, " and parameters method are: ",
            paste0(names(parameters.method), " = ", parameters.method, collapse = ", "))
    cat("\n")
  } else {
    message("Downsampling method is: ", d.method)
    cat("\n")
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
             trans <- flowCore::estimateLogicle(fcs, channels = flowCore::colnames(fcs), m = 5.5)
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

    } else if (d.method == "density") {
      all.exprs <- downsamplingDensity(file = file,
                                       fcs = fcs,
                                       exprs.raw = exprs.raw,
                                       exprs.sub = exprs.sub,
                                       exclude.pctile = parameters.method$exclude.pctile,
                                       target.pctile = parameters.method$target.pctile,
                                       target.percent = parameters.method$target.percent,
                                       target.number = parameters.method$target.number)

      exprs.raw <- all.exprs$exprs.raw
      exprs.sub <- all.exprs$exprs.sub
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
    samples.r <- c(samples.r, rep(sample, nrow(exprs.raw)))
    samples <- c(samples, rep(sample, nrow(exprs.sub)))

  }

  if (!is.null(exclude.markers)) {
    exprs.r <- exprs.r[, !colnames(exprs.r) %in% exclude.markers]
    exprs   <- exprs[, !colnames(exprs) %in% exclude.markers]
  }

  res <- methods::new("UMAPdata",
                      matrix.expression.r = exprs.r,
                      matrix.expression = exprs,
                      samples = samples,
                      raw.markers = colnames(exprs),
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

# @title Internal - Computes the downsampling density-based
#
# @description This function aims to perform a downsampling density-based cells
#
# @param file a character vector specifying the path of the tab-separated or FCS files to load
# @param fcs
# @param exprs.raw a data.frame providing the raw marker expressions
# @param exprs.sub a data.frame providing the marker expressions
# @param exclude.pctile a numeric value specifying the density threshold to be excluded
# @param target.pctile a numeric value specifying the density threshold to be maintained
# @param target.percent a numeric value providing the percentage of cells to downsample for each sample
# @param target.number a numeric value providing the number of cells to downsample for each sample
#
# @return a list containing the two dowsampled expression data.frame
#
# @export
#
downsamplingDensity <- function(file,
                                fcs,
                                exprs.raw,
                                exprs.sub,
                                exclude.pctile,
                                target.pctile,
                                target.percent,
                                target.number) {

  idxs <- names(flowCore::markernames(fcs))
  density <- spade::SPADE.density(exprs.sub[, idxs])

  ### from spade
  exprs.sub <- cbind(exprs.sub, "density" = density)

  d.idxs <- match("density", colnames(exprs.sub))
  boundary <- stats::quantile(exprs.sub[, d.idxs],
                              c(exclude.pctile, target.pctile), names = FALSE)

  idx.e <- exprs.sub[, d.idxs] > boundary[1]
  exprs.sub <- exprs.sub[idx.e, ]
  exprs.raw <- exprs.raw[idx.e, ]

  if (!is.null(target.percent) && is.null(target.number)) {
    target.number <- round(target.percent * nrow(exprs.sub))
    message("Number of cells: ", target.number, " for ", basename(file))
  }

  density <- exprs.sub[, d.idxs]
  if (is.null(target.number) && is.null(target.percent)) {
    boundary <- boundary[2]
    idx.c <- boundary / density > stats::runif(nrow(exprs.sub))
    exprs.sub <- exprs.sub[idx.c, ]
    exprs.raw <- exprs.raw[idx.c, ]

  } else if (target.number < nrow(exprs.sub)) {
    density.s <- sort(density)
    cdf <- rev(cumsum(1 / rev(density.s)))
    boundary <- target.number / cdf[1]
    if (boundary > density.s[1]) {
      targets <- (target.number - seq(1, length(density.s))) / cdf
      boundary <- targets[which.min(targets - density.s > 0)]
    }

    idx.c <- boundary / density > stats::runif(length(density))
    exprs.sub <- exprs.sub[idx.c, ]
    exprs.raw <- exprs.raw[idx.c, ]
  } else if (target.number > nrow(exprs.sub)) {
    stop("Number of events is too high")
  }
  ### from spade

  exprs.sub <- subset(exprs.sub, select = -c(get("density")))

  return(list("exprs.raw" = exprs.raw, "exprs.sub" = exprs.sub))
}

#' @title Renames markers within a UMAPdata object
#'
#' @description This function aims to rename cell markers stored within a UMAPdata object.
#'
#' This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
#'
#' @param UMAPdata a UMAPdata object
#' @param marker.names a character vector providing the new marker names to use
#'
#' @return a S4 object of class 'UMAPdata'
#'
#'@export
#'
renameMarkers <- function(UMAPdata,
                         marker.names) {

  if (length(marker.names) != length(unique(marker.names))) {
    stop("The same marker is used twice")
  }

  if (ncol(UMAPdata@matrix.expression) != length(marker.names)) {
    stop("The number of markers given is different from the initial number of markers")
  }

  colnames(UMAPdata@matrix.expression) <- marker.names
  validObject(UMAPdata)

  return(UMAPdata)
}

#' @title Assigns meta-information about biological samples
#'
#' @description This function aims to attach meta-information to each biological sample.
#'
#' Especially, the biological individual, the biological condition and the time point of each sample can be specified for subsequent analyses.
#'
#' @param UMAPdata a UMAPdata object
#' @param metadata a data.frame containing contextual information about the biological samples. This data.frame must have 3 columns specifying for each sample the associated individual (column named 'individual'), the biological condition (column named 'condition') and the time point (column named 'timepoint')
#'
#' @return a S4 object of class 'UMAPdata'
#'
#' @export
#'
assignMetadata <- function(UMAPdata,
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

  if (all(unique(UMAPdata@samples) != rownames(metadata))) {
    stop("sample names (rownames) of the metadata are inconcistent with the sample names stored in the UMAPdata object")
  }

  UMAPdata@metadata <- metadata

  validObject(UMAPdata)
  return(UMAPdata)
}
