# Avoid R CMD check notes for pipe operator
utils::globalVariables(c("."))

# Avoid R CMD check notes for variables used in dplyr/tidyr operations
utils::globalVariables(c(
  "MGI.symbol", "HGNC.symbol", 
  "Chromosome.scaffold.name", "Chromosome.scaffold.name.1",
  "mouse_gene_symbol", "human_gene_symbol", 
  "mouse_chromosome", "human_chromosome",
  "Gene.Name", "samples", "exprs"
))

#' Helper function to get biomart dictionary
#' 
#' Internal function to retrieve the biomart dictionary from the package environment,
#' loading it if necessary
#' 
#' @return A data frame containing gene mappings
#' @keywords internal
get_biomart_dict <- function() {
  # Check if biomart_dict exists in package environment
  if (!exists("biomart_dict", envir = .strpip_env)) {
    # Try to load from local file
    tryCatch({
      import_biomart_local()
    }, error = function(e) {
      stop("BioMart dictionary not available. ", 
           "Run import_biomart() or import_biomart_local() first.", 
           call. = FALSE)
    })
  }
  
  # Return the dictionary
  get("biomart_dict", envir = .strpip_env)
}

#' Validate organism parameter
#' 
#' Internal function to validate the organism parameter
#' 
#' @param org Organism string to validate
#' @return A validated organism string ("human" or "mouse")
#' @keywords internal
validate_organism <- function(org) {
  # Convert to lowercase
  org <- tolower(org)
  
  # Validate
  if (!org %in% c("human", "mouse")) {
    stop("Invalid organism. Must be 'human' or 'mouse'", call. = FALSE)
  }
  
  return(org)
} 