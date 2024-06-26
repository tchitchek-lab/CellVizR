#' @title Computes differential analysis statistics for cell clusters
#'
#' @description This function aims to identify of differentially abundant clusters.
#'
#' Such clusters correspond to cell clusters having abundances statistically different between two biological conditions.
#' The statistical test used for the comparisons can be defined by users.
#' For each cluster, the p-value, log2 fold-change and effect size relative to the reference condition are computed.
#' Statistical comparison can be performed in a paired and unpaired manner.
#' Computed p-values can be corrected for multiple testing.
#'
#' @param Celldata a Celldata object
#' @param condition a character value providing the name of the condition to be compared
#' @param ref.condition a character value providing the name of reference condition
#' @param comparison a character value providing the name of comparison
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcoxon' or 't-test'
#' @param p.adjust a character value providing the type of p-value adjustment method to use. Possible values are: 'none', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param paired a boolean value indicating if individuals are paired
#'
#' @return a S4 object of class 'Celldata'
#'
#' @export
#' @import rstatix
#'
computeStatistics <- function(Celldata,
                              condition,
                              ref.condition,
                              comparison = "cmp",
                              test.statistics = c("wilcox.test", "t.test"),
                              p.adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                              paired = FALSE) {

  test.statistics <- match.arg(test.statistics)
  p.adjust <- match.arg(p.adjust)

  checkmate::qassert(condition, "S+")
  checkmate::qassert(ref.condition, "S+")
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(p.adjust, "S1")
  checkmate::qassert(paired, "B1")

  message("Computing of the ", test.statistics, " for: ", paste(condition,collapse = ","), " vs. ", paste(ref.condition,collapse = ","))

  if (comparison %in% unique(Celldata@statistic$comparison)) {
    stop("The statistic slot already exists")
  }

  cmp_effsize <- function(test.statistics, paired, x, y) {
    data <- rbind(cbind(x, "x"), cbind(y, "y"))
    data <- data.frame(data)
    colnames(data) <- c("value", "grp")
    data$value <- as.numeric(data$value)

    test.eff <- ifelse(test.statistics == "t.test", "rstatix::cohens_d", "rstatix::wilcox_effsize")
    effsize <- do.call(eval(parse(text = test.eff)), list(data = data,
                                                          stats::as.formula("value~grp"),
                                                          ci = FALSE, paired = paired))
    statistic <- effsize$effsize
    statistic <- as.numeric(statistic)

    return(statistic)
  }

  abundances <- Celldata@matrix.abundance
  abundances <- data.matrix(abundances)

  stats <- data.frame()
  for (cluster in as.character(rownames(abundances))) {
    values.condition <- abundances[cluster, condition, drop = TRUE]
    values.ref.condition <- abundances[cluster, ref.condition, drop = TRUE]

    pv <- do.call(test.statistics, list(x = values.condition,
                                        y = values.ref.condition,
                                        exact = FALSE, paired = paired))$p.value
    effsize <- cmp_effsize(test.statistics, paired = paired,
                           x = values.condition, y = values.ref.condition)

    lfc <- log(mean(values.condition) / mean(values.ref.condition)) / log(2)
    effsize <- effsize * sign(lfc)

    stats <- rbind(stats, cbind(cluster = cluster, pv = pv,
                               lfc = lfc, effsize = effsize))
  }

  stats <- data.frame(stats)
  colnames(stats) <- c("clusters", "pvalue", "lfc", "effsize")
  stats$pvalue <- as.numeric(stats$pvalue)
  stats$lfc <- as.numeric(stats$lfc)
  stats$effsize <- as.numeric(stats$effsize)
  
  if(p.adjust != "none"){
	stats$pvalue = p.adjust(stats$pvalue,method = p.adjust)
  }
  
  if (nrow(Celldata@statistic) == 0) {
    Celldata@statistic <- cbind(comparison, stats)
  } else {
    Celldata@statistic <- rbind(Celldata@statistic, cbind(comparison, stats))
  }

  validObject(Celldata)
  return(Celldata)
  
}
