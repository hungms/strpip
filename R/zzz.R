#' @importFrom utils packageVersion
#'
#' @title Internal Package Environment
#' @description This environment is used to store package-level data that needs to persist
#' between function calls without polluting the user's global environment.
#' @keywords internal
.strpip_env <- new.env(parent = emptyenv())

#' Package startup message
#'
#' Displays the package version when the package is attached
#' @param libname The library where the package is installed
#' @param pkgname The name of the package
#' @importFrom utils packageVersion
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("strpip", packageVersion("strpip")))
}

#' Initialize package when loaded
#' 
#' This function is called when the package is loaded via library() or require()
#' @param libname The library where the package is installed
#' @param pkgname The name of the package
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL
  
  packages <- c(
    "dplyr", "tidyr", "tibble", "stringr", "utils", "rlang","magrittr")

  invisible(lapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
    }
  }))

  if (requireNamespace("OmnipathR", quietly = TRUE)) {
    library("OmnipathR")}
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    library("biomaRt")}
}
