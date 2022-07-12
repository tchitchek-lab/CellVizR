#' @title Performs the upsampling of downsampled events
#'
#' @description This function aims to perform upsample downsampled events based on an existing UMAPdata object and existing cell events stored in tab-separated or FCS files. 
#' 
#' Importantly, the identification of cell clusters must have been performed prior to this operation. 
#'  
#' @param UMAPdata a UMAPdata object
#' @param files a character vector providing the path of the tab-separated or FCS files
#'  
#' @return a S4 object of class 'UMAPdata'
performUpsampling <- function(UMAPdata, 
                              files) {
  
  checkmate::qassert(files, "S*")
  
  sample_files <- basename(files)
  sample_files <- gsub(".fcs","",sample_files)
  sample_files <- gsub(".FCS","",sample_files)
  sample_files <- gsub(".txt","",sample_files)
  sample_files <- gsub(".TXT","",sample_files)
  
  files = files[sample_files %in% unique(UMAPdata@samples)]
  UMAPdata_wodownsampling <- import(files)
  
  downsampled.exp <- UMAPdata@matrix.expression
  colnames(downsampled.exp) <- UMAPdata@raw.markers
  
  full.exp <- UMAPdata_wodownsampling@matrix.expression
  full.exp <- full.exp[,colnames(full.exp) %in% colnames(downsampled.exp)]
  
  chk.downsampled <- apply(downsampled.exp,1,sum)
  chk.upsampled  <- apply(full.exp,1,sum)
  
  upsampled.exp <- full.exp[!chk.upsampled %in% chk.downsampled,]
  upsampled.samples <- UMAPdata_wodownsampling@samples[!chk.upsampled %in% chk.downsampled]
  upsampled.clusters <- UMAPdata_wodownsampling@identify.clusters[!chk.upsampled %in% chk.downsampled]
  
  downsampled.exp$cluster <- as.numeric(UMAPdata@identify.clusters) ######
  downsampled.centers <- plyr::ddply(downsampled.exp,"cluster",function(x){
    x$cluster <- NULL
    centers <- apply(x,2,stats::median,rm.na=TRUE)
    return(centers)
  })
  downsampled.centers <- downsampled.centers[order(downsampled.centers$cluster),]
  downsampled.centers$cluster <- NULL
  downsampled.exp$cluster <- NULL
  
  knn <- FNN::knnx.index(downsampled.centers, upsampled.exp, k=1, algorithm="kd_tree")
  
  UMAPdata@matrix.expression  <- rbind(downsampled.exp, upsampled.exp)
  UMAPdata@samples            <- c(UMAPdata@samples, upsampled.samples)
  UMAPdata@identify.clusters  <- c(UMAPdata@identify.clusters, knn)
  
  tofill <- data.frame(dim1 = rep(NA,nrow(upsampled.exp)), dim2 = rep(NA,nrow(upsampled.exp)))
  UMAPdata@manifold <- rbind(UMAPdata@manifold, tofill)
  samples <- UMAPdata@samples
  
  message(paste0("computing cell cluster count matrix..."))
  UMAPdata@matrix.cell.count <- computeCellCounts(samples <- UMAPdata@samples,
                                                  proj <- UMAPdata@matrix.expression,
                                                  clusters <- UMAPdata@identify.clusters)
  
  message(paste0("computing cell cluster abundance matrix..."))
  UMAPdata@matrix.abundance <- computeClusterAbundances(count <- UMAPdata@matrix.cell.count)
  
  return(UMAPdata)
  
}
