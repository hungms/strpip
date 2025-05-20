#' Import BioMart data from Ensembl archive
#'
#' Retrieves archived ensembl biomarts from web and creates a mapping dictionary between
#' human and mouse genes.
#'
#' @param host URL to retrieve archived ensembl biomarts. Default is the December 2021 archive.
#' @return A data frame containing gene mappings and chromosome information
#' @examples
#' \dontrun{
#' biodict_mouse <- import_biomart_mouse()
#' head(biodict_mouse)
#' }
#' @export
import_biomart_mouse <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {

   if(local){
      biodict_mouse <- read.table(file.path(system.file("extdata", package = "strpip"), paste0("biodict_mouse_release-", release, ".tsv")), header = TRUE, sep = "\t")
      assign("biodict_mouse", biodict_mouse, envir = .strpip_env)
      get("biodict_mouse", envir = .strpip_env)}
   else{
      mouse_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl")
      attributes <- c("mgi_symbol", "chromosome_name")
      biodict_mouse <- getBM(attributes = attributes, mart = mouse_biomart)
      colnames(biodict_mouse) <- c("mouse_gene_symbol", "mouse_chromosome")
      return(biodict_mouse)
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
#' biodict_mouse <- import_biomart_mouse()
#' head(biodict_mouse)
#' }
#' @export

import_biomart_human <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {

   if(local){
      biodict_human <- read.table(file.path(system.file("extdata", package = "strpip"), paste0("biodict_human_release-", release, ".tsv")), header = TRUE, sep = "\t")
      assign("biodict_human", biodict_human, envir = .strpip_env)
      get("biodict_human", envir = .strpip_env)}
   else{
      human_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl")
      attributes <- c("hgnc_symbol", "chromosome_name")
      biodict_human <- getBM(attributes = attributes, mart = human_biomart)
      colnames(biodict_human) <- c("human_gene_symbol", "human_chromosome")
      return(biodict_human)}
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
#' biodict <- import_biomart()
#' head(biodict)
#' }
#' @export
import_biomart <- function(host = 'https://dec2021.archive.ensembl.org', local = TRUE, release = "105") {
   if(local) {
      # Get the package's extdata directory
      pkg_dir <- system.file("extdata", package = "strpip")
      if (pkg_dir == "") {
         stop("Package directory not found. Is strpip installed correctly?")
      }
      
      # Construct file path
      file_path <- file.path(pkg_dir, paste0("biodict_release-", release, ".tsv"))
      if (!file.exists(file_path)) {
         stop("Biomart dictionary file not found at: ", file_path)
      }
      
      # Read the file with explicit encoding and error handling
      biodict <- tryCatch({
         utils::read.table(file_path, 
                         header = TRUE, 
                         sep = "\t",
                         fileEncoding = "UTF-8",
                         stringsAsFactors = FALSE)
      }, error = function(e) {
         stop("Error reading biomart dictionary file: ", e$message, "\nFile path: ", file_path)
      })
      
      # Validate the data
      if (nrow(biodict) == 0) {
         stop("Biomart dictionary file is empty")
      }
      
      # Store in package environment and return
      assign("biodict", biodict, envir = .strpip_env)
      get("biodict", envir = .strpip_env)
   } else {
      # Create biomart connections
      human_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl")
      mouse_biomart <- biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl")
      
      # Get linked datasets
      biodict_mouse <- biomaRt::getLDS(attributes = c("mgi_symbol", "chromosome_name"), # Mouse attributes
                     filters = "", # No specific filter, retrieve all genes
                     values = NULL, # No specific values, retrieve all genes
                     mart = mouse_biomart, # Mouse BioMart
                     attributesL = c("hgnc_symbol", "chromosome_name"), # Human attributes
                     martL = human_biomart, # Human BioMart
                     uniqueRows = F) # Remove duplicate rows
      biodict_human <- biomaRt::getLDS(attributes = c("hgnc_symbol", "chromosome_name"), # Mouse attributes
                     filters = "", # No specific filter, retrieve all genes
                     values = NULL, # No specific values, retrieve all genes
                     mart = human_biomart, # Mouse BioMart
                     attributesL = c("mgi_symbol", "chromosome_name"), # Human attributes
                     martL = mouse_biomart, # Human BioMart
                     uniqueRows = F) # Remove duplicate rows
      
      # Process the data
      biodict_mouse <- biodict_mouse %>% 
         dplyr::mutate(
            mouse_gene_symbol = MGI.symbol,
            human_gene_symbol = HGNC.symbol,
            mouse_chromosome = Chromosome.scaffold.name,
            human_chromosome = Chromosome.scaffold.name.1) %>%
         dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
         dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") | stringr::str_detect(human_gene_symbol, "[A-Z]"))

      # Process the data
      biodict_human <- biodict_human %>% 
         dplyr::mutate(
            human_gene_symbol = HGNC.symbol,
            mouse_gene_symbol = MGI.symbol,
            human_chromosome = Chromosome.scaffold.name,
            mouse_chromosome = Chromosome.scaffold.name.1) %>%
         dplyr::select(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome) %>%
         dplyr::filter(stringr::str_detect(mouse_gene_symbol, "[a-z]") | stringr::str_detect(human_gene_symbol, "[A-Z]"))
      
      biodict <- rbind(biodict_mouse, biodict_human) %>%
            distinct(mouse_gene_symbol, mouse_chromosome, human_gene_symbol, human_chromosome, .keep_all = TRUE)
      return(biodict)
   }
}