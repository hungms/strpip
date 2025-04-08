#' Get genes from X and Y chromosomes
#'
#' Retrieves all gene symbols located on the X and Y chromosomes for a specified organism.
#'
#' @param org Organism to query. Either "human" or "mouse". Default is "human".
#' @param ... Additional arguments passed to `import_biomart_human` or `import_biomart_mouse`.
#' @return A character vector containing gene symbols from X and Y chromosomes.
#' @examples
#' \dontrun{
#' xy_genes <- get_xy_genes("human")
#' head(xy_genes)
#' }
#' @export
get_xy_genes <- function(org = "human", ...) {
  # Validate organism parameter
  org <- validate_organism(org)
  
  # Get the biomart dictionary
  if(org == "human"){
    biomart_dict_human <- import_biomart_human(...)}
  else{
    biomart_dict_mouse <- import_biomart_mouse(...)}
  
  # Extract genes on X and Y chromosomes
  xy_genes <- get(paste0("biomart_dict_", org)) %>% 
      dplyr::filter(!!rlang::sym(paste0(org, "_chromosome")) %in% c("X", "Y")) %>%
      dplyr::select(!!rlang::sym(paste0(org, "_gene_symbol"))) %>%
      .[[1]] %>%
      unique()
      
  # Check if we found any genes
  if (length(xy_genes) == 0) {
    warning("No X/Y chromosome genes found for organism: ", org, call. = FALSE)
  }
    
  return(xy_genes)
}

#' Get genes from mitochondrial chromosome
#'
#' Retrieves all gene symbols located on the mitochondrial chromosome for a specified organism.
#'
#' @param org Organism to query. Either "human" or "mouse". Default is "human".
#' @param ... Additional arguments passed to `import_biomart_human` or `import_biomart_mouse`.
#' @return A character vector containing mitochondrial gene symbols.
#' @examples
#' \dontrun{
#' mt_genes <- get_mt_genes("human")
#' head(mt_genes)
#' }
#' @export
get_mt_genes <- function(org = "human", ...) {
  # Validate organism parameter
  org <- validate_organism(org)
  
  # Get the biomart dictionary
  if(org == "human"){
    biomart_dict_human <- import_biomart_human(...)}
  else{
    biomart_dict_mouse <- import_biomart_mouse(...)}
  
  # Extract mitochondrial genes
  mt_genes <- get(paste0("biomart_dict_", org)) %>% 
      dplyr::filter(!!rlang::sym(paste0(org, "_chromosome")) %in% c("MT")) %>%
      dplyr::select(!!rlang::sym(paste0(org, "_gene_symbol"))) %>%
      .[[1]] %>%
      unique()
      
  # Check if we found any genes
  if (length(mt_genes) == 0) {
    warning("No mitochondrial genes found for organism: ", org, call. = FALSE)
  }
    
  return(mt_genes)
}

#' Get genes matching specific patterns
#'
#' Retrieves gene symbols matching predefined patterns for specific gene families
#' (BCR, TCR, MHC, HB, RB, MT) for a specified organism.
#'
#' @param org Organism to query. Either "human" or "mouse". Default is "human".
#' @param str Character vector of patterns to search for. Valid options are:
#'   "bcr", "tcr", "mhc", "hb", "rb", "mt". Default is "rb".
#' @return A character vector containing matching gene symbols.
#' @examples
#' \dontrun{
#' rb_genes <- get_str_genes("human", "rb")
#' head(rb_genes)
#' }
#' @export
get_str_genes <- function(org = "human", str = c("rb"), ...) {
  # Validate organism parameter
  org <- validate_organism(org)
  
  # Validate pattern parameter
  valid_patterns <- c("bcr", "tcr", "mhc", "hb", "rb", "mt")
  if (!all(str %in% valid_patterns)) {
    invalid_patterns <- str[!str %in% valid_patterns]
    stop("Invalid pattern(s): ", paste(invalid_patterns, collapse = ", "), 
         ". Valid patterns are: ", paste(valid_patterns, collapse = ", "), 
         call. = FALSE)
  }
  
  # Combine patterns
  pattern_strings <- paste0(unlist(
    mget(paste0(str, ".string"), 
         envir = asNamespace("strpip"), 
         inherits = FALSE)), 
    collapse = "|")
  
  # Get the biomart dictionary
  if(org == "human"){
    biomart_dict_human <- import_biomart_human(...)}
  else{
    biomart_dict_mouse <- import_biomart_mouse(...)}
  
  # Extract genes matching the pattern
  genes <- get(paste0("biomart_dict_", org)) %>% 
      dplyr::filter(stringr::str_detect(
        !!rlang::sym(paste0(org, "_gene_symbol")), 
        pattern_strings)) %>%
      dplyr::select(!!rlang::sym(paste0(org, "_gene_symbol"))) %>%
      .[[1]] %>%
      unique()
      
  # Check if we found any genes
  if (length(genes) == 0) {
    warning("No genes found matching patterns: ", paste(str, collapse = ", "), 
            " for organism: ", org, call. = FALSE)
  }
    
  return(genes)
}

#' Regular expression pattern for BCR genes
#'
#' Pattern to identify B-cell receptor genes
#' @export
bcr.string <- "^I[Gg][HKLhkl][VDJCAEMGLvdjcaemgl]|^AC233755"

#' Regular expression pattern for TCR genes
#'
#' Pattern to identify T-cell receptor genes
#' @export
tcr.string <- "^T[Rr][ABCDGabcdg][VDJCvdjc]"

#' Regular expression pattern for hemoglobin genes
#'
#' Pattern to identify hemoglobin genes
#' @export
hb.string <- "^H[B][ABDEGMPQZ]?\\d*$|^H[b][abdegmpqz]?\\d*"

#' Regular expression pattern for MHC genes
#'
#' Pattern to identify major histocompatibility complex genes
#' @export
mhc.string <- "^HLA-|^H2-"

#' Regular expression pattern for ribosomal genes
#'
#' Pattern to identify ribosomal protein genes
#' @export
rb.string <- "^R[Pp][SsLl]"

#' Regular expression pattern for mitochondrial genes
#'
#' Pattern to identify mitochondrial genes
#' @export
mt.string <- "^[Mm][Tt]-"