#' Add Gene Annotations
#'
#' Adds functional annotations to a data frame containing gene symbols using pre-computed
#' OmniPath annotations. The function supports both human and mouse gene symbols.
#'
#' @param df Data frame containing gene symbols
#' @param gene_column Name of the column containing gene symbols. Default is "gene"
#' @param org.from Organism type. Either "human" or "mouse". Default is "human"
#' @param org.to Target organism. Either "human" or "mouse". Default is "mouse"
#' @param release Ensembl release version to use for annotations. Default is "105"
#' @param one.to.many Logical. Default is FALSE to return only unique mappings (one/many to one). If TRUE, returns all possible mouse gene mappings (one to many), including cases where one human gene maps to multiple mouse genes.
#' @return A data frame with additional columns containing functional annotations:
#'         - Pathway information
#'         - Protein-protein interactions
#'         - Gene regulatory relationships
#'         - Other functional annotations from OmniPath
#' @importFrom utils read.table
#' @importFrom dplyr %>% left_join filter select mutate
#' @export
run_annotation <- function(df, gene_column = "gene", org.from = "human", org.to = "mouse", 
                         release = "105", one.to.many = FALSE) {

    # Validate inputs and ensure data frame format
    df <- validate_df(df, gene_column)
    validate_org(org.from, genes = df[[gene_column]])
    validate_org(org.to)

    # Construct annotation file path
    annotation_file <- system.file(
        "extdata", 
        paste0("omnipathdb_", org.from, "_release-", release, ".tsv"),
        package = "strpip")

    # Read annotation data
    annot <- read.table(annotation_file, sep = "\t", header = TRUE)
    
    # Rename the gene column in annotations to match the desired join column
    names(annot)[names(annot) == paste0(org.from, "_gene_symbol")] <- gene_column
    
    # Merge annotations with input data
    result <- df %>%
        dplyr::left_join(annot, by = gene_column)
    
    # If one.to.many is FALSE, remove multiple mappings
    if (!one.to.many) {
        target_col <- paste0(org.to, "_gene_symbol")
        result[[target_col]] <- gsub(",.*", "", result[[target_col]])}
    
    return(result)
}