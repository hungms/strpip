#' Get BioMart Dictionary
#'
#' Internal function to access or load the biomart dictionary.
#' This ensures the dictionary is available without polluting the global environment.
#'
#' @return A data frame containing gene mappings between human and mouse
#' @keywords internal
get_biomart_dict <- function() {
  if (!exists("biomart_dict", envir = .strpip_env)) {
    import_biomart_local()
  }
  get("biomart_dict", envir = .strpip_env)
}

#' Check if Gene Symbols Are Valid
#'
#' Internal function to check if gene symbols are valid for a given organism.
#'
#' @param genes Character vector of gene symbols to check
#' @param org Organism type, either "human" or "mouse"
#' @return Logical indicating whether the genes are valid
#' @keywords internal
validate_gene_symbols <- function(genes, org = "human") {
  if (!is.character(genes)) {
    stop("Gene symbols must be provided as a character vector", call. = FALSE)
  }
  
  if (org == "mouse" && !any(grepl("[a-z]", genes))) {
    warning("Unexpected gene format for mouse: mouse genes typically contain lowercase letters", call. = FALSE)
    return(FALSE)
  }
  
  if (org == "human" && !any(grepl("[A-Z]", genes))) {
    warning("Unexpected gene format for human: human genes typically contain uppercase letters", call. = FALSE)
    return(FALSE)
  }
  
  return(TRUE)
}

#' Validate Organism Parameter
#'
#' Internal function to validate the organism parameter.
#'
#' @param org Organism to validate
#' @return The validated organism or an error if invalid
#' @keywords internal
validate_organism <- function(org) {
  if (!org %in% c("human", "mouse")) {
    stop("Organism must be either 'human' or 'mouse'", call. = FALSE)
  }
  return(org)
} 