% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{renameMarkers}
\alias{renameMarkers}
\title{Renames markers within a Celldata object}
\usage{
renameMarkers(Celldata, marker.names)
}
\arguments{
\item{Celldata}{a Celldata object}

\item{marker.names}{a character vector providing the new marker names to use}
}
\value{
a S4 object of class 'Celldata'
}
\description{
This function aims to rename cell markers stored within a Celldata object.

This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
}
