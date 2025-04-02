#' strpip: Data Structure Pipelines for Genomic Analysis
#'
#' @description
#' A comprehensive R package for genomic data structure pipelines, providing tools for:
#' \itemize{
#'   \item Data format conversions (GCT, GMT, TSV)
#'   \item Gene symbol mapping between species
#'   \item Gene set operations and filtering
#'   \item Gene expression data processing
#'   \item Functional annotation integration
#'   \item Specialized gene set identification (XY, MT, BCR, TCR, etc.)
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{Data Format Operations}{`read_gmt()`, `read_gct()`, `write_gmt()`, `write_gct()`}
#'   \item{Gene Symbol Mapping}{`convert_mouse_to_human()`, `convert_human_to_mouse()`}
#'   \item{Gene Set Operations}{`get_xy_genes()`, `get_mt_genes()`, `get_str_genes()`}
#'   \item{Expression Processing}{`summarize_genes()`}
#'   \item{Annotation}{`run_annotation()`}
#' }
#'
#' @section Dependencies:
#' The package depends on:
#' \itemize{
#'   \item biomaRt: For accessing Ensembl BioMart
#'   \item dplyr, tidyr, tibble, stringr: For data manipulation
#'   \item magrittr: For pipe operations
#'   \item OmniPathR: For pathway data
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr tidyr tibble stringr
#' @importFrom magrittr %>%
#' @importFrom utils packageVersion data
#' @importFrom rlang .data
## usethis namespace: end
NULL
