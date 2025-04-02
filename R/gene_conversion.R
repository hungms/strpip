#' Import BioMart data from Ensembl archive
#'
#' Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
#' human and mouse genes.
#'
#' @param host URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.
#' @return A data frame containing gene mappings and chromosome information
#' @examples
#' \dontrun{
#' biomart_dict <- import_biomart()
#' head(biomart_dict)
#' }
#' @export
import_biomart <- function(host = 'https://dec2021.archive.ensembl.org') {
   # Create biomart connections
   tryCatch({
     human_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl")
     mouse_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl")
     
     # Get linked datasets
     biomart_dict <- biomaRt::getLDS(attributes = c("mgi_symbol", "chromosome_name"), # Mouse attributes
                    filters = "", # No specific filter, retrieve all genes
                    values = NULL, # No specific values, retrieve all genes
                    mart = mouse_biomart, # Mouse BioMart
                    attributesL = c("hgnc_symbol", "chromosome_name"), # Human attributes
                    martL = human_biomart, # Human BioMart
                    uniqueRows = TRUE) # Remove duplicate rows
     
     # Process the data
     biomart_dict <- biomart_dict %>% 
       dplyr::mutate(
           mouse_gene_symbol = MGI.symbol,
           human_gene_symbol = HGNC.symbol,
           mouse_chromosome = Chromosome.scaffold.name,
           human_chromosome = Chromosome.scaffold.name.1) %>%
       dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
       dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") | stringr::str_detect(human_gene_symbol, "[A-Z]"))
     
     # Store in package environment
     assign("biomart_dict", biomart_dict, envir = .strpip_env)
     
     # Return for direct use
     return(biomart_dict)
   }, error = function(e) {
     stop("Error accessing BioMart: ", e$message, 
          "\nEnsure you have a working internet connection or try using import_biomart_local() instead.", 
          call. = FALSE)
   })
}

#' Import BioMart data from local file
#'
#' Retrieves gene mapping data from a local TSV file instead of querying Ensembl.
#' This is faster than `import_biomart()` but requires the data file to be present.
#'
#' @param release Ensembl release number. Default is "105".
#' @return A data frame containing gene mappings between species
#' @examples
#' \dontrun{
#' biomart_dict <- import_biomart_local()
#' head(biomart_dict)
#' }
#' @export
import_biomart_local <- function(release = "105") {
    file <- paste0(system.file("extdata", package = "strpip"), "/release-", release, ".tsv")
    
    if (!file.exists(file)) {
      stop("BioMart data file not found: ", file, 
           "\nMake sure the package is properly installed.", call. = FALSE)
    }
    
    tryCatch({
      biomart_dict <- utils::read.table(file, header = TRUE, sep = "\t") %>%
        dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") & 
                     stringr::str_detect(human_gene_symbol, "[A-Z]"))
      
      # Store in package environment
      assign("biomart_dict", biomart_dict, envir = .strpip_env)
      
      # Return for direct use
      return(biomart_dict)
    }, error = function(e) {
      stop("Error reading BioMart data file: ", e$message, call. = FALSE)
    })
}

#' Convert mouse gene symbols to human gene symbols
#' 
#' Converts a list of mouse (MGI) gene symbols to their human (HGNC) equivalents.
#' By default, returns only one-to-one mappings to ensure accuracy.
#'
#' @param genes A vector of mouse gene symbols to convert
#' @param one.to.many Logical. If TRUE, returns all possible human gene mappings,
#'   including cases where one mouse gene maps to multiple human genes.
#' @return If `one.to.many = FALSE`, returns a vector of unique human gene symbols.
#'   If `one.to.many = TRUE`, returns a data frame with all possible mappings.
#' @examples
#' \dontrun{
#' mouse_genes <- c("Trp53", "Cd4", "Cd8a")
#' human_genes <- convert_mouse_to_human(mouse_genes)
#' print(human_genes)
#' }
#' @export
convert_mouse_to_human <- function(genes, one.to.many = FALSE) {
   # Validate input
   if (!is.character(genes)) {
     stop("Gene list must be a character vector", call. = FALSE)
   }
   
   # Get the biomart dictionary
   biomart_dict <- get_biomart_dict()
   
   # Filter to relevant genes and create mapping
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::filter(mouse_gene_symbol %in% genes) %>%
      dplyr::select(c(human_gene_symbol, mouse_gene_symbol)) %>% 
      dplyr::distinct(human_gene_symbol, .keep_all = TRUE)
      
   # Return appropriate format based on parameters
   if (one.to.many) {
      return(biomart_dict_tmp)
   } else {
      human.genes <- biomart_dict_tmp$human_gene_symbol %>% unique()
      return(human.genes)
   }
}

#' Convert human gene symbols to mouse gene symbols
#' 
#' Converts a list of human (HGNC) gene symbols to their mouse (MGI) equivalents.
#' By default, returns only one-to-one mappings to ensure accuracy.
#'
#' @param genes A vector of human gene symbols to convert
#' @param one.to.many Logical. If TRUE, returns all possible mouse gene mappings,
#'   including cases where one human gene maps to multiple mouse genes.
#' @return If `one.to.many = FALSE`, returns a vector of unique mouse gene symbols.
#'   If `one.to.many = TRUE`, returns a data frame with all possible mappings.
#' @examples
#' \dontrun{
#' human_genes <- c("TP53", "CD4", "CD8A")
#' mouse_genes <- convert_human_to_mouse(human_genes)
#' print(mouse_genes)
#' }
#' @export
convert_human_to_mouse <- function(genes, one.to.many = FALSE) {
   # Validate input
   if (!is.character(genes)) {
     stop("Gene list must be a character vector", call. = FALSE)
   }
   
   # Get the biomart dictionary
   biomart_dict <- get_biomart_dict()
   
   # Filter to relevant genes and create mapping
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::filter(human_gene_symbol %in% genes) %>%
      dplyr::select(c(human_gene_symbol, mouse_gene_symbol)) %>% 
      dplyr::distinct(mouse_gene_symbol, .keep_all = TRUE)
      
   # Return appropriate format based on parameters
   if (one.to.many) {
      return(biomart_dict_tmp)
   } else {
      mouse.genes <- biomart_dict_tmp$mouse_gene_symbol %>% unique()
      return(mouse.genes)
   }
}

