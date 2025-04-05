
#' Import BioMart data from Ensembl archive
#'
#' Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
#' human and mouse genes.
#'
#' @param host URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.
#' @return A data frame containing gene mappings and chromosome information
#' @examples
#' \dontrun{
#' biomart_dict_mouse <- import_biomart_mouse()
#' head(biomart_dict_mouse)
#' }
#' @export
import_biomart_mouse <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {

   if(local){
      biomart_dict_mouse <- read.table(file.path(system.file("extdata", package = "strpip"), paste0("release-", release, "_biomart_dict_mouse.tsv")), header = TRUE, sep = "\t")
      assign("biomart_dict_mouse", biomart_dict_mouse, envir = .strpip_env)
      get("biomart_dict_mouse", envir = .strpip_env)}
   else{
      mouse_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl")
      attributes <- c("mgi_symbol", "chromosome_name")
      biomart_dict_mouse <- getBM(attributes = attributes, mart = mouse_biomart)
      colnames(biomart_dict_mouse) <- c("mouse_gene_symbol", "mouse_chromosome")
      return(biomart_dict_mouse)
      }
}

#' Import BioMart data from Ensembl archive
#'
#' Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
#' human and mouse genes.
#'
#' @param host URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.
#' @return A data frame containing gene mappings and chromosome information
#' @examples
#' \dontrun{
#' biomart_dict_mouse <- import_biomart_mouse()
#' head(biomart_dict_mouse)
#' }
#' @export

import_biomart_human <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {

   if(local){
      biomart_dict_human <- read.table(file.path(system.file("extdata", package = "strpip"), paste0("release-", release, "_biomart_dict_human.tsv")), header = TRUE, sep = "\t")
      assign("biomart_dict_human", biomart_dict_human, envir = .strpip_env)
      get("biomart_dict_human", envir = .strpip_env)}
   else{
      human_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl")
      attributes <- c("hgnc_symbol", "chromosome_name")
      biomart_dict_human <- getBM(attributes = attributes, mart = human_biomart)
      colnames(biomart_dict_human) <- c("human_gene_symbol", "human_chromosome")
      return(biomart_dict_human)}
}

#' Import BioMart data from Ensembl archive
#'
#' Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
#' human and mouse genes.
#'
#' @param host URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.
#' @return A data frame containing gene mappings and chromosome information
#' @examples
#' \dontrun{
#' biomart_dict <- import_biomart_orthologs()
#' head(biomart_dict)
#' }
#' @export
import_biomart_orthologs <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {

   if(local){
      biomart_dict_orthologs <- utils::read.table(file.path(system.file("extdata", package = "strpip"), paste0("/release-", release, "_biomart_dict_orthologs.tsv")), header = TRUE, sep = "\t")
      assign("biomart_dict_orthologs", biomart_dict_orthologs, envir = .strpip_env)
      get("biomart_dict_orthologs", envir = .strpip_env)}
   else{
      # Create biomart connections
      tryCatch({
         human_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl")
         mouse_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl")
         
         # Get linked datasets
         biomart_dict_mouse_to_human <- biomaRt::getLDS(attributes = c("mgi_symbol", "chromosome_name"), # Mouse attributes
                        filters = "", # No specific filter, retrieve all genes
                        values = NULL, # No specific values, retrieve all genes
                        mart = mouse_biomart, # Mouse BioMart
                        attributesL = c("hgnc_symbol", "chromosome_name"), # Human attributes
                        martL = human_biomart, # Human BioMart
                        uniqueRows = F) # Remove duplicate rows
         biomart_dict_human_to_mouse <- biomaRt::getLDS(attributes = c("hgnc_symbol", "chromosome_name"), # Mouse attributes
                        filters = "", # No specific filter, retrieve all genes
                        values = NULL, # No specific values, retrieve all genes
                        mart = human_biomart, # Mouse BioMart
                        attributesL = c("mgi_symbol", "chromosome_name"), # Human attributes
                        martL = mouse_biomart, # Human BioMart
                        uniqueRows = F) # Remove duplicate rows
         
         # Process the data
         biomart_dict_mouse_to_human <- biomart_dict_mouse_to_human %>% 
            dplyr::mutate(
               mouse_gene_symbol = MGI.symbol,
               human_gene_symbol = HGNC.symbol,
               mouse_chromosome = Chromosome.scaffold.name,
               human_chromosome = Chromosome.scaffold.name.1) %>%
            dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
            dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") | stringr::str_detect(human_gene_symbol, "[A-Z]"))

         # Process the data
         biomart_dict_human_to_mouse <- biomart_dict_human_to_mouse %>% 
            dplyr::mutate(
               human_gene_symbol = HGNC.symbol,
               mouse_gene_symbol = MGI.symbol,
               human_chromosome = Chromosome.scaffold.name,
               mouse_chromosome = Chromosome.scaffold.name.1) %>%
            dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
            dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") | stringr::str_detect(human_gene_symbol, "[A-Z]"))
         
         biomart_dict_orthologs <- rbind(biomart_dict_mouse_to_human, biomart_dict_human_to_mouse) %>%
               distinct(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome, .keep_all = TRUE)
         return(biomart_dict_orthologs)
            
         }, error = function(e) {
         stop("Error accessing BioMart: ", e$message, 
               "\nEnsure you have a working internet connection or try using import_biomart_local() instead.", 
               call. = FALSE)
         })}
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