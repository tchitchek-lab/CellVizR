#' @title Computes differential analysis statistics for cell clusters 
#'
#' @description This function aims to compute the statistics of Differentially Abundant Clusters 'DAC'. 
#' 
#' DAC correspond to cell clusters having abundances statistically different between two biological conditions.
#' The statistical test used for the comparisons can be defined by users. 
#' For each cluster, the p-value, log2 fold-change and effect size relative to the reference condition are computed.  
#' Statistical comparison can be performed in a paired and unpaired manner. 
#'  
#' @param UMAPdata a UMAPdata object  
#' @param condition a character value providing the name of the condition to be compared 
#' @param ref.condition a character value providing the name of reference condition
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcoxon' or 't-test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#'  
#' @return a S4 object of class 'UMAPdata'
#' 
#' @export 
#' 
computeStatistics = function(UMAPdata, 
                             condition, 
                             ref.condition,
                             test.statistics = c("wilcoxon", "t-test"),
                             paired = c("paired", "unpaired")) {
  
  message(paste0("Computing: ", condition, " vs. ", ref.condition))
  
  abundances = UMAPdata@matrix.abundance
  abundances = data.matrix(abundances)
  
  clusters = rownames(abundances)
  
  stats = data.frame()
  for(cluster in as.character(clusters)) {
    values.condition = abundances[cluster, grepl(condition, colnames(abundances)), drop = TRUE]
    values.ref.condition = abundances[cluster, grepl(ref.condition, colnames(abundances)), drop = TRUE]
    if(length(values.condition)<=2 || length(values.ref.condition)<=2) {
      pv = 1
    }else{
      pv = stats::wilcox.test(values.condition, values.ref.condition, exact = FALSE)$p.value
    }
    
    if(all(values.condition==0) || all(values.ref.condition==0)) {
      statistic = 0
    }else if(length(values.condition)<=2 || length(values.ref.condition)<=2) {
      statistic = 0
    }else{
      forES = rbind(cbind(values.condition, "A"), cbind(values.ref.condition, "B"))
      forES = data.frame(forES)
      colnames(forES) = c("value", "grp")
      forES$value = as.numeric(forES$value)
      
      tryCatch({
        effsize = rstatix::wilcox_effsize(data = forES, value~grp, ci=FALSE)
        statistic = effsize$effsize
        statistic = as.numeric(statistic)
      }, error = function(e) {
        statistic = 0
      })
    }
    
    lfc = log(mean(values.condition)/mean(values.ref.condition))/log(2)
    statistic = statistic*sign(lfc)
    
    stats = rbind(stats, cbind(cluster=cluster, pv=pv, lfc=lfc, statistic=statistic))
  }
  
  stats = data.frame(stats)
  colnames(stats) = c("clusters", "pvalue", "lfc", "statistic")
  stats$pvalue = as.numeric(stats$pvalue)
  stats$lfc = as.numeric(stats$lfc)
  stats$statistic = as.numeric(stats$statistic)
  
  patients = colnames(abundances)
  patients = patients[grepl(paste0(condition, "|", ref.condition), patients)]
  #patients = as.vector(sapply(patients, function(x) paste(strsplit(x, "_")[[1]][c(2)], collapse = "_")))
  stats$n.patients = length(unique(patients))
  
  return(stats)
  
}


