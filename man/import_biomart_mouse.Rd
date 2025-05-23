% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart_database.R
\name{import_biomart_mouse}
\alias{import_biomart_mouse}
\title{Import BioMart data from Ensembl archive}
\usage{
import_biomart_mouse(
  host = "https://dec2021.archive.ensembl.org",
  local = TRUE,
  release = "105"
)
}
\arguments{
\item{host}{URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.}
}
\value{
A data frame containing gene mappings and chromosome information
}
\description{
Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
human and mouse genes.
}
\examples{
\dontrun{
biodict_mouse <- import_biomart_mouse()
head(biodict_mouse)
}
}
