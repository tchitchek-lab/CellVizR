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
                             test.statistics = c("wilcox.test", "t.test"),
                             paired = FALSE) {
  
  test.statistics <- match.arg(test.statistics)
  
  checkmate::qassert(condition, "S1")
  checkmate::qassert(ref.condition, "S+")
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(paired, "B1")
  
  message(paste0("Computing: ", condition, " vs. ", ref.condition))
  
  comparison = paste0(condition," vs. ",ref.condition)
  if(comparison %in% unique(UMAPdata@statistic$comparison)){
	stop("!!!")
  }
	
  cmp_effsize <- function(test.statistics, paired, x, y){	
		data           <- rbind(cbind(x, "x"), cbind(y, "y"))
		data           <- data.frame(data)
		colnames(data) <- c("value", "grp")
		data$value     <- as.numeric(data$value)
		
		test.eff       <- ifelse(test.statistics=="t.test","rstatix::cohens_d","rstatix::wilcox_effsize")
		effsize        <- do.call(eval(parse(text=test.eff)),list(data = data, stats::as.formula("value~grp"), ci=FALSE, paired=paired))
		statistic      <- effsize$effsize
		statistic      <- as.numeric(statistic)
		
		return(statistic)
	} 
	
  abundances = UMAPdata@matrix.abundance
  abundances = data.matrix(abundances)
  
  stats = data.frame()
  for(cluster in as.character(rownames(abundances))) {
  
    values.condition     = abundances[cluster, grepl(condition, colnames(abundances)), drop = TRUE]
    values.ref.condition = abundances[cluster, grepl(ref.condition, colnames(abundances)), drop = TRUE]
	
	pv = do.call(test.statistics, list(x=values.condition, y=values.ref.condition, exact = FALSE, paired=paired))$p.value
    effsize = cmp_effsize(test.statistics, paired=paired, x=values.condition, y=values.ref.condition)
    
	lfc = log(mean(values.condition)/mean(values.ref.condition))/log(2)
    effsize = effsize*sign(lfc)
    
    stats = rbind(stats, cbind(cluster=cluster, pv=pv, lfc=lfc, effsize=effsize))
  }
  
  stats = data.frame(stats)
  colnames(stats) = c("clusters", "pvalue", "lfc", "effsize")
  stats$pvalue    = as.numeric(stats$pvalue)
  stats$lfc       = as.numeric(stats$lfc)
  stats$effsize   = as.numeric(stats$effsize)
  
  if(nrow(UMAPdata@statistic)==0){
	UMAPdata@statistic = cbind(comparison,stats)
  }else{
	UMAPdata@statistic = rbind(UMAPdata@statistic,cbind(comparison,stats))
  }
  
  return(UMAPdata)
  
}


