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
#' @param filetype a character vector specifying xxx
#' @param transform a character value containing the type of the transformation to apply. Possible values are: 'logicle', 'arcsinh', 'logarithmic' or 'none' 
#' @param downsampling a numeric value providing the number of cells to downsample for each sample
#' @param d.method a character value containing the type of the dowsampling to apply. Possible values are: 'uniform' or 'density'
#' @param exclude.markers a character vector providing the marker names to be excluded during the import
#' @param seed a numeric value providing the random seed to use during stochastic operations 
#' 
#' @return a S4 object of class 'UMAPdata'
#' 
#' @export
import = function(files, 
                  filetype="fcs",
                  transform = c("logicle","arcsinh", "logarithmic", "none"), 
                  downsampling = NULL,
                  d.method = c("uniform","density"), 
                  exclude.markers = NULL,
                  seed = 42) {
  
  transform <- match.arg(transform)
  d.method <- match.arg(d.method)
  
  checkmate::qassert(filetype, "S1")
  checkmate::qassert(transform, "S1")
  checkmate::qassert(downsampling, c("0", "N1"))
  # checkmate::qassert(d.method, "S1")
  checkmate::qassert(exclude.markers, c("0", "S*"))
  checkmate::qassert(seed, "N1")
  
  exprs <- data.frame()
  samples <- c()
  
  for(file in files) {
    
    message(paste0("importing ",basename(file)," file"))
    
    if(filetype == "fcs"){
      fcs <- flowCore::read.FCS(file)
    }else{
      exprs <- utils::read.delim(file,sep="\t")
      exprs <- exprs[,!colnames(exprs) %in% exclude.markers]
      fcs <- suppressWarnings(createFlowframe(exprs))
    }
    
    switch(transform, 
           logicle = {
             trans <- flowCore::estimateLogicle(fcs, channels = flowCore::colnames(fcs), m = 5.5)
             fcs <- flowCore::transform(fcs, trans)
           },
           arcsinh = {
             trans.arcsinh <- flowCore::arcsinhTransform(a=0, b=0.2)
             marker.trans <- flowCore::colnames(fcs)
             trans <- flowCore::transformList(marker.trans, trans.arcsinh)
             fcs <- flowCore::transform(fcs, trans)
           },
           logarithmic = {
             trans.log <- flowCore::logTransform(logbase = 10, r=1,d=1)
             marker.trans <- flowCore::colnames(fcs)
             trans <- flowCore::transformList(marker.trans, trans.log)
             fcs <- flowCore::transform(fcs, trans)
           }, 
           none = {
             fcs <- fcs 
           })
    
    exprs.sub <- flowCore::exprs(fcs)
    exprs.sub <- exprs.sub[,colnames(exprs.sub) %in% names(flowCore::markernames(fcs))]
    
    if(!is.null(downsampling)){
      set.seed(seed)
      exprs.sub <- exprs.sub[sample(nrow(exprs.sub),min(downsampling,nrow(exprs.sub))),]
    }
    
    colnames(exprs.sub) <- flowCore::markernames(fcs)
    
    exprs <- rbind(exprs, exprs.sub)
    
    sample <- gsub(".fcs", "", basename(file))
    samples <- c(samples, rep(sample, nrow(exprs.sub)))
    
  }
  
  if(!is.null(exclude.markers)){
    exprs <- exprs[,!colnames(exprs) %in% exclude.markers]
  }
  
  res <- methods::new("UMAPdata",
                      matrix.expression = exprs,
                      samples = samples,
                      raw.markers = colnames(exprs),
                      matrix.abundance = data.frame())
  
  return(res)
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
#' @export
#' 
renameMarkers = function(UMAPdata,
                         marker.names){
  
  checkmate::qassert(marker.names, "S*")
  
  # Ne permet pas d'afficher la bonne erreur 
  if(all(marker.names != unique(marker.names))){
    stop("xxx")
  }
  if(length(colnames(UMAPdata@matrix.expression)) != length(marker.names)) {
    stop("xxx")
  }
  
  colnames(UMAPdata@matrix.expression) <- marker.names
  
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
assignMetadata = function(UMAPdata, 
                          metadata) {
  
  checkmate::qassert(metadata, c("D2", "D3"))
  
  if(length(unique(colnames(metadata)))!=length(colnames(metadata))){
    stop("colnames are not unique")
  }
  
  if(!all(colnames(metadata) %in% c("individual", "condition", "timepoint"))) {
    stop("colnames must be 'individual', 'condition', or 'timepoint'")
  }
  
  if(!all(c("individual") %in% colnames(metadata))) {
    stop("colnames 'individual' is missing")
  }
  
  if(all(unique(UMAPdata@samples) != rownames(metadata))){
    stop("sample names (rownames) of the metadata are inconcistent with the sample names stored in the UMAPdata object")
  }
  
  UMAPdata@metadata <- metadata
  
  return(UMAPdata)
}
