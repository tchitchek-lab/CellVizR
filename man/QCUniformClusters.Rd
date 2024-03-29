% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{QCUniformClusters}
\alias{QCUniformClusters}
\title{Computes the percentage of clusters with uniform phenotypes}
\usage{
QCUniformClusters(
  Celldata,
  uniform.test = c("both", "uniform", "IQR"),
  th.pvalue = 0.05,
  th.IQR = 2,
  plot.device = TRUE
)
}
\arguments{
\item{Celldata}{a Celldata object}

\item{uniform.test}{a character providing the name of test assessment to perform. Possible value are : 'both', 'uniform', 'IQR'}

\item{th.pvalue}{a numeric value providing the p-value threshold of the Hartigan's dip test (unimodal if pvalue > th.pvalue)}

\item{th.IQR}{a numeric value providing the IQR (interquartile range) threshold to assume a distribution as uniform}

\item{plot.device}{a boolean value specifying if result representation must be displayed}
}
\value{
a numerical value corresponding to the percentage of cell cluster with unimodal expression and a low spread
}
\description{
This function aims to identify and show cell clusters having a uniform phenotype.

A uniform cluster corresponds to a cluster that have a unimodal expression and a low spread of expression for all its markers.
}
\details{
-'uniform' corresponds to the verification of the unimodal distribution of markers with a Hartigans test

-'IQR' corresponds to the verification of the distribution of markers so that they are not below the IQR threshold (interquantile range)

-'both' corresponds to the combination of the two parameters : uniform and IQR
}
