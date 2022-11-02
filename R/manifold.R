#' @title Generates a manifold of cell events
#'
#' @description This function aims to generate a manifold representation for cell events stored in a Celldata object
#'
#' This function allows the use of several non-linear dimension reduction techniques such as UMAP, t-SNE or LargeVis
#' The whole set of cell markers or specific cell markers can be used during the dimensionality reduction process
#'
#' @param Celldata a Celldata object
#' @param type a character value specifying the type of manifold to compute. Possible values are: 'UMAP' for Uniform Manifold Approximation and Projection, 'tSNE' for t Stochastic Neighbor Embedding, and 'lvish' for LargeVis
#' @param markers a character vector providing the cell markers to use for the manifold generation
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param verbose a boolean value indicating if computational details must be displayed on the console
#' @param ... Other arguments passed on to methods
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#'
generateManifold <- function(Celldata,
                             type = c("UMAP", "tSNE", "lvish"),
                             markers = NULL,
                             seed = 42,
                             verbose = TRUE,
                             ...) {

  type <- match.arg(type)

  checkmate::qassert(type, "S1")
  checkmate::qassert(markers, c("0", "S*"))
  checkmate::qassert(seed, "N1")
  checkmate::qassert(verbose, "B1")

  matrix.exprs <- Celldata@matrix.expression

  if (is.null(markers)) {
    markers <- names(matrix.exprs)
  }

  matrix.exprs.subset <- matrix.exprs[, names(matrix.exprs)
                                      %in% markers]

  message("Manifold markers are: ", paste0(unlist(markers), collapse = ", "))
  message("Manifold method is: ", type)
  cat("\n")

  switch(type,
         UMAP = {
           generate.manifold <- generateManifoldUMAP(matrix.exprs.subset,
                                                     seed = seed,
                                                     verbose = verbose,
                                                     ...)
         },
         tSNE = {
           generate.manifold <- generateManifoldtSNE(matrix.exprs.subset,
                                                     seed = seed,
                                                     verbose = verbose,
                                                     ...)
         },
         lvish = {
           generate.manifold <- generateManifoldlvish(matrix.exprs.subset,
                                                      seed = seed,
                                                      verbose = verbose,
                                                      ...)
         })

  Celldata@manifold <- generate.manifold

  manifold.params <- list(type = type,
                          markers = markers,
                          verbose = verbose,
                          ...)
  Celldata@manifold.params <- manifold.params

  validObject(Celldata)
  return(Celldata)
}

# @title Internal - Generates a UMAP manifold of cell events
#
# @description This function is used internally a manifold representation based on the UMAP algorithms
#
# @param exprs a data.frame providing the marker expressions
# @param n.neighbors a numeric value providing the size of local neighborhood used for manifold approximation (please refer to the function 'umap' of the 'uwot' package)
# @param n.components a numeric value providing the dimension of the space to embed into (please refer to the function 'umap' of the 'uwot' package)
# @param metric a character value providing type of distance metric to use to find nearest neighbors (please refer to the function 'umap' of the 'uwot' package)
# @param n.epochs a numeric value providing number of epochs to use during the optimization of the embedded coordinates (please refer to the function 'umap' of the 'uwot' package)
# @param n.threads a numeric value providing the number of threads to use (please refer to the function 'umap' of the 'uwot' package)
# @param n.sgd.threads a numeric value providing the number of threads to use during stochastic gradient descent (please refer to the function 'umap' of the 'uwot' package)
# @param scale a character value providing the scaling to apply (please refer to the function 'umap' of the 'uwot' package)
# @param seed a numeric value providing the random seed to use in stochastic operation
# @param verbose a boolean value indicating if computational details must be displayed on the console
# @param ... Other arguments passed on to methods
#
# @return a data.frame containing the manifold coordinates
#
generateManifoldUMAP <- function(exprs,
                                 seed,
                                 verbose,
                                 ...) {

  do.call("set.seed", list(seed))
  umap <- uwot::umap(exprs,
                     verbose = verbose,
                     ...)

  umap <- data.frame(umap)
  colnames(umap) <- c("dim1", "dim2")

  return(umap)
}

# @title Internal - Generates a t-SNE manifold of cell events
#
# @description This function is used internally to compute a manifold representation based on the t-SNE algorithms
#
# @param exprs a data.frame containing the marker expressions
# @param dims a numeric value providing the output dimensionality (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param initial.dims a numeric value providing the number of dimensions that should be retained in the initial PCA step (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param perplexity a numeric value providing perplexity parameter (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param theta a numeric value providing speed/accuracy trade-off (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param max.iter a numeric value providing number of iterations (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param seed a numeric value providing the random seed to use in stochastic operation
# @param verbose a boolean value indicating if computational details must be displayed on the console (please refer to the function 'Rtsne' of the 'Rtsne' package)
# @param ... Other arguments passed on to methods
#
# @return a data.frame containing the manifold coordinates
#
generateManifoldtSNE <- function(exprs,
                                 seed,
                                 verbose,
                                 ...) {

  do.call("set.seed", list(seed))
  tSNE <- Rtsne::Rtsne(exprs,
                       verbose = verbose,
                       ...)

  tSNE <- tSNE$Y
  tSNE <- data.frame(tSNE)
  colnames(tSNE) <- c("dim1", "dim2")

  return(tSNE)
}

# @title Internal - Generates a LargeVis manifold of cell events
#
# @description This function is used internally to compute a manifold representation based on the LargeVis algorithms
#
# @param exprs a data.frame containing the marker expressions
# @param perplexity a numeric value providing the size of the local neighborhood used (please refer to the function 'lvish' of the 'uwot' package)
# @param n.neighbors a numeric value providing the number of neighbors to use when calculating the perplexity (please refer to the function 'lvish' of the 'uwot' package)
# @param n.components a numeric value providing the dimension of the space to embed into (please refer to the function 'lvish' of the 'uwot' package)
# @param metric a character value providing type of distance metric to use to find nearest neighbors (please refer to the function 'lvish' of the 'uwot' package)
# @param n.epochs a numeric value providing number of epochs to use during the optimization of the embedded coordinates (please refer to the function 'lvish' of the 'uwot' package)
# @param n.threads a numeric value providing the number of threads to use (please refer to the function 'lvish' of the 'uwot' package)
# @param n.sgd.threads a numeric value providing the number of threads to use during stochastic gradient descent (please refer to the function 'lvish' of the 'uwot' package)
# @param scale a character value providing the scaling to apply (please refer to the function 'lvish' of the 'uwot' package)
# @param seed a numeric value providing the random seed to use in stochastic operation
# @param verbose a boolean value indicating if computational details must be displayed on the console
# @param ... Other arguments passed on to methods
#
# @return a data.frame containing the manifold coordinates
#
generateManifoldlvish <- function(exprs,
                                  seed,
                                  verbose,
                                  ...) {

  do.call("set.seed", list(seed))
  lvish <- uwot::lvish(exprs,
                       verbose = verbose,
                       ...)

  lvish <- data.frame(lvish)
  colnames(lvish) <- c("dim1", "dim2")

  return(lvish)
}
