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
#' Converts a vector of mouse (MGI) gene symbols to their human (HGNC) equivalents.
#' By default, returns a data frame with all possible mappings.
#'
#' @param genes A vector of mouse gene symbols to convert
#' @return A data frame with all possible mappings.
#' @examples
#' \dontrun{
#' mouse_genes <- c("Trp53", "Cd4", "Cd8a")
#' human_genes <- convert_mouse_to_human(mouse_genes)
#' print(human_genes)
#' }
#' @export
convert_mouse_to_human <- function(genes) {

   # Get the biomart dictionary
   biomart_dict <- get_biomart_dict()
   
   # Filter to relevant genes and create mapping
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::filter(mouse_gene_symbol %in% genes) %>%
      dplyr::select(c(mouse_gene_symbol, human_gene_symbol))
      
   # Return appropriate format based on parameters
   return(biomart_dict_tmp)
}

#' Convert human gene symbols to mouse gene symbols
#' 
#' Converts a list of human (HGNC) gene symbols to their mouse (MGI) equivalents.
#' By default, returns only one-to-one mappings to ensure accuracy.
#'
#' @param genes A vector of human gene symbols to convert
#' @param one.to.many Logical. Default is FALSE to return only unique mappings (one/many to one). If TRUE, returns all possible mouse gene mappings (one to many), including cases where one human gene maps to multiple mouse genes.
#' @return If `one.to.many = FALSE`, returns a vector of unique mouse gene symbols.
#'   If `one.to.many = TRUE`, returns a data frame with all possible mappings.
#' @examples
#' \dontrun{
#' human_genes <- c("TP53", "CD4", "CD8A")
#' mouse_genes <- convert_human_to_mouse(human_genes)
#' print(mouse_genes)
#' }
#' @export
convert_human_to_mouse <- function(genes){
   
   # Get the biomart dictionary
   biomart_dict <- get_biomart_dict()
   
   # Filter to relevant genes and create mapping
   biomart_dict_tmp <- biomart_dict %>% 
      dplyr::filter(human_gene_symbol %in% genes) %>%
      dplyr::select(c(human_gene_symbol, mouse_gene_symbol))

   return(biomart_dict_tmp)}



#' Convert orthologs between human and mouse for a vector of genes
#'
#' Converts a vector of gene symbols between human and mouse species.
#'
#' @param genes A vector of gene symbols to convert
#' @param mode The direction of conversion, either "human_to_mouse" or "mouse_to_human"
#' @param one.to.many Logical. Default is FALSE to return only unique mappings (one/many to one). If TRUE, returns all possible mouse gene mappings (one to many), including cases where one human gene maps to multiple mouse genes.
#' @return A vector of converted gene symbols
#' @examples
#' \dontrun{
#' mouse_genes <- c("Trp53", "Cd4", "Cd8a")
#' human_genes <- convert_orthologs_vector(mouse_genes, mode = "mouse_to_human", one.to.many = TRUE)
#' print(human_genes)
#' }
#' @export
convert_orthologs_vector <- function(genes, mode = "human_to_mouse", one.to.many = TRUE){
   stopifnot(mode %in% c("human_to_mouse", "mouse_to_human"))
   stopifnot(is.logical(one.to.many))

   if(mode == "human_to_mouse"){
      orthologs_df <- convert_human_to_mouse(genes)}
   else{
      orthologs_df <- convert_mouse_to_human(genes)}

   if(one.to.many){
      mapped_genes <- orthologs_df %>% .[[2]] %>% unique(.)}
   else{
      mapped_genes <- orthologs_df %>% distinct(.[[1]], .keep_all = TRUE) %>% .[[2]] %>% unique(.)}
   return(mapped_genes)
}

#' Convert orthologs between human and mouse for an gene expression matrix
#'
#' Converts a matrix of gene symbols between human and mouse species.
#'
#' @param matrix A matrix of gene symbols to convert
#' @param mode The direction of conversion, either "human_to_mouse" or "mouse_to_human"
#' @param many.to.one Logical. Default is FALSE to return only unique mappings (one/many to one). If TRUE, returns all possible mouse gene mappings (one to many), including cases where one human gene maps to multiple mouse genes.
#' @return A matrix of converted gene symbols
#' @examples
#' \dontrun{
#' matrix <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' rownames(matrix) <- c("gene1", "gene2", "gene3", "gene4", "gene5")
#' colnames(matrix) <- c("sample1", "sample2", "sample3", "sample4", "sample5")
#'
#' human_matrix <- convert_orthologs_matrix(matrix, mode = "human_to_mouse", many.to.one = TRUE)
#' mouse_matrix <- convert_orthologs_matrix(matrix, mode = "mouse_to_human", many.to.one = TRUE)
#' }
#' @export
convert_orthologs_matrix <- function(matrix, mode = "human_to_mouse", many.to.one = TRUE, normalized = FALSE){

   stopifnot(mode %in% c("human_to_mouse", "mouse_to_human"))
   stopifnot(is.logical(many.to.one))
   stopifnot(any(stringr::str_detect(rownames(matrix), "[A-Z]")))

   if(mode == "human_to_mouse"){
      orthologs_df <- convert_human_to_mouse(rownames(matrix))}
   else{
      orthologs_df <- convert_mouse_to_human(rownames(matrix))}

   orig.species <- colnames(orthologs_df)[[1]]
   mapped.species <- colnames(orthologs_df)[[2]]

   if(many.to.one){
      mapped_matrix <- matrix %>%
         as.data.frame() %>%
         rownames_to_column("gene") %>%
         pivot_longer(!gene, names_to = "sample", values_to = "exprs") %>%
         merge(., orthologs_df %>% distinct(!!sym(orig.species), .keep_all = TRUE), by.x = "gene", by.y = orig.species, all.x = TRUE) %>%
         filter(!!sym(mapped.species) != "NA") %>%
         pivot_wider(names_from = sample, values_from = exprs) %>%
         dplyr::select(-c(gene)) %>%
         as.data.frame(.)

      mapped_vector <- mapped_matrix[[mapped.species]]
      mapped_matrix[[mapped.species]] <- NULL
      mapped_matrix <- summarize_genes(mapped_matrix, gene_sym_vec = mapped_vector, normalized = normalized)
         }
   else{
      orthologs_df <- orthologs_df %>% 
         group_by(human_gene_symbol) %>%
         mutate(hn = n()) %>%
         group_by(mouse_gene_symbol) %>%
         mutate(mn = n()) %>%
         filter(hn == 1 & mn == 1)
      orthologs_df <- orthologs_df[c(1,2)]
      mapped_matrix <- matrix[which(rownames(matrix) %in% orthologs_df[[orig.species]]),]

      mapped_matrix <- mapped_matrix %>%
         as.data.frame() %>%
         merge(., orthologs_df, by.x = 0, by.y = orig.species, all.x = TRUE) %>%
         dplyr::select(-c("Row.names")) %>%
         column_to_rownames(mapped.species) %>%
         as.data.frame(.)
         }
   return(mapped_matrix)
}

