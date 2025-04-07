#' Add Gene Annotations
#'
#' Adds functional annotations to a data frame containing gene symbols using pre-computed
#' OmniPath annotations. The function supports both human and mouse gene symbols.
#'
#' @param df Data frame containing gene symbols
#' @param gene_column Name of the column containing gene symbols. Default is "gene"
#' @param org Organism type. Either "human" or "mouse". Default is "human"
#' @param release Ensembl release version to use for annotations. Default is "105"
#' @return A data frame with additional columns containing functional annotations:
#'         - Pathway information
#'         - Protein-protein interactions
#'         - Gene regulatory relationships
#'         - Other functional annotations from OmniPath
#' @examples
#' \donttest{
#' # Create example data frame with human genes
#' expression_df <- data.frame(
#'   gene = c("TP53", "BRCA1", "EGFR"),
#'   expression = c(1.2, 0.8, 2.1)
#' )
#' 
#' # Annotate human genes
#' annotated_df <- run_annotation(expression_df, gene_column = "gene", org = "human")
#' 
#' # Create example data frame with mouse genes
#' mouse_df <- data.frame(
#'   gene = c("Trp53", "Brca1", "Egfr"),
#'   expression = c(1.2, 0.8, 2.1)
#' )
#' 
#' # Annotate mouse genes
#' annotated_df <- run_annotation(mouse_df, gene_column = "gene", org = "mouse")
#' }
#' @export
run_annotation <- function(df, gene_column = "gene", org = "human", release = "105") {
    # Check if OmnipathR is available
    if (!requireNamespace("OmnipathR", quietly = TRUE)) {
        warning("The OmnipathR package is needed for annotation functions. ",
                "Please install it with: BiocManager::install('OmnipathR')",
                call. = FALSE)
        return(df)
    }
    
    # Validate organism
    org <- validate_organism(org)
    
    # Validate data frame
    if (!is.data.frame(df)) {
        stop("Input must be a data frame", call. = FALSE)
    }
    
    # Validate gene column
    if (!gene_column %in% colnames(df)) {
        stop("Column '", gene_column, "' not found in data frame", call. = FALSE)
    }
    
    # Validate gene symbols match organism type
    genes <- df[[gene_column]]
    if (org == "mouse" && !any(grepl("[a-z]", genes))) {
        warning("Mouse gene symbols typically contain lowercase letters. Check if your gene symbols are correct.", call. = FALSE)
    } else if (org == "human" && !any(grepl("[A-Z]", genes))) {
        warning("Human gene symbols typically contain uppercase letters. Check if your gene symbols are correct.", call. = FALSE)
    }
    
    # Find annotation file
    files <- list.files(system.file("extdata", package = "strpip"), full.names = TRUE)
    pattern <- paste0("release-", release, "_", org, "_omnipath_db.tsv")
    annotation_file <- files[grep(pattern, files)][1]
    
    # Read and merge annotations
    tryCatch({
        # Read annotation data
        annot <- utils::read.table(annotation_file, header = TRUE, sep = "\t", 
                                 stringsAsFactors = FALSE)
        
        # Merge with input data frame
        result_df <- merge(df, annot, 
                         by.x = gene_column, 
                         by.y = paste0(org, "_gene_symbol"), 
                         all.x = TRUE)
        
        return(result_df)
    }, error = function(e) {
        warning("Error while annotating genes: ", e$message, 
               ". Returning original data frame.", call. = FALSE)
        return(df)
    })
}