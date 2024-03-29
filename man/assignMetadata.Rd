% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{assignMetadata}
\alias{assignMetadata}
\title{Assigns meta-information to biological samples}
\usage{
assignMetadata(Celldata, metadata)
}
\arguments{
\item{Celldata}{a Celldata object}

\item{metadata}{a data.frame containing contextual information about the biological samples. This data.frame must have 3 columns specifying for each sample the associated individual (column named 'individual'), the biological condition (column named 'condition') and the time point (column named 'timepoint'). Rownames must correspond to biological samples imported within the Celldata object.}
}
\value{
a S4 object of class 'Celldata'
}
\description{
This function aims to attach meta-information to biological samples.

For each biological sample, the biological individual, the biological condition and the time point can be specified for subsequent analyses.
}
