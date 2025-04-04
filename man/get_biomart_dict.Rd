% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-biomart.R, R/utils.R
\name{get_biomart_dict}
\alias{get_biomart_dict}
\title{Get BioMart Dictionary}
\usage{
get_biomart_dict()

get_biomart_dict()
}
\value{
A data frame containing gene mappings between human and mouse

A data frame containing gene mappings
}
\description{
Internal function to access or load the biomart dictionary.
This ensures the dictionary is available without polluting the global environment.

Internal function to retrieve the biomart dictionary from the package environment,
loading it if necessary
}
\keyword{internal}
