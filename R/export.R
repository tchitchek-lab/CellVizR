# @title Internal - Creates of a FlowFrame object
#
# @description This function is used internally to create a flowframe object with the purpose extracting marker intensities to FCS files. 
#
# @param intensities a data.frame providing the cell profile intensities
#
# @return a flowframe object containing the marker expression
# 
createFlowframe <- function(intensities) {
  
  markers <- colnames(intensities)
  p <- c()    
  description <- list() 
  
  description[["$DATATYPE"]] <- "F"
  
  for (i in seq(1,ncol(intensities))) {
    
    name  <- markers[i]
    min   <- min(intensities[,i])
    max   <- max(intensities[,i])
    range <- max-min+1
    
    l           <- matrix(c(name,name,range,min,max),nrow=1)
    colnames(l) <- c("name","desc","range","minRange","maxRange")
    rownames(l) <- paste0("$P",i) 
    p           <- rbind(p,l)
    
    description[[paste("$P",i,"N",sep="")]] <- name;
    description[[paste("$P",i,"S",sep="")]] <- name;
    description[[paste("$P",i,"R",sep="")]] <- toString(range);
    description[[paste("$P",i,"B",sep="")]] <- "32";
    description[[paste("$P",i,"E",sep="")]] <- "0,0";
  }
  
  intensities <- as.matrix(intensities)
  dataframe <- methods::as(data.frame(p), "AnnotatedDataFrame")
  
  flowframe <- flowCore::flowFrame(intensities, dataframe, description=description)
  
  return(flowframe)
}

#' @title Exports cell expression profiles to TSV or FCS files 
#'
#' @description Exports cell expresion profiles from a UMAPdata object to a tab-separated or FCS files.
#' 
#' Cell expression profiles can be exported for a set of given samples and for a set of given cell clusters. 
#'  
#' @param UMAPdata a UMAPdata object 
#' @param filename a character value providing the name of the output file
#' @param clusters a character vector containing the identifiers of the cell clusters to export. By default, all clusters are extracted
#' @param samples a character vector containing the names of biological samples to export. By default, all samples are extracted
#' 
#' @return none
#' 
#' @name export
#' @rdname export-methods
methods::setGeneric("export",function(UMAPdata,
                                      filename,
                                      clusters = NULL,
                                      samples = NULL) { standardGeneric("export") })



#' @rdname export-methods
#' @export
methods::setMethod("export",c("UMAPdata"),
                   function(UMAPdata,filename,clusters=NULL,samples = NULL){
                     
                     checkmate::qassert(filename, c("0","S1"))
                     checkmate::qassert(clusters, c("0","S*"))
                     checkmate::qassert(samples, c("0","S*"))
                     
                     exprs <- UMAPdata@matrix.expression
                     
                     if(length(UMAPdata@manifold)!=0){
                       exprs <- cbind(exprs,UMAPdata@manifold)
                     }
                     if(length(UMAPdata@identify.clusters)!=0){
                       exprs <- cbind(exprs,cluster=as.numeric(UMAPdata@identify.clusters))
                     }
                     
                     exprs <- cbind(exprs,samples=UMAPdata@samples)
                     
                     
                     if(!is.null(samples)){
                       exprs <- exprs[exprs$samples %in% samples,]
                     }
                     if(!is.null(clusters)){
                       exprs <- exprs[exprs$clusters %in% clusters,]
                     }
                     
                     exprs <- data.frame(exprs)
                     
                     if(any(tools::file_ext(filename)%in% c("fcs","FCS"))){
                       
                       exprs$samples <- as.numeric(factor(exprs$samples))					 
                       
                       flowFrame <- suppressWarnings(createFlowframe(exprs))
                       flowCore::write.FCS(flowFrame,filename=filename)
                       
                     }else{
                       
                       exprs$sample.id <- as.numeric(factor(exprs$samples))	
                       
                       utils::write.table(exprs, file = filename, sep="\t", row.names=FALSE)
                     } 
                     
                   }
)

