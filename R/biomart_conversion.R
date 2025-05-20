#' Convert mouse gene symbols
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
get_conversion_dict <- function(org.from, org.to) {
    # Validate inputs
    validate_org(org.from)
    validate_org(org.to)
    stopifnot(org.from != org.to)
    
    # Get the biomart dictionary
    biodict <- import_biomart()
    
    # Get column names
    source_col <- paste0(org.from, "_gene_symbol")
    target_col <- paste0(org.to, "_gene_symbol")
    
    # Filter and create mapping
    biodict <- biodict %>%
        dplyr::filter(!is.na(!!sym(target_col))) %>%
        dplyr::filter(grepl("[A-Z]", !!sym(target_col))) %>%
        dplyr::select(!!sym(source_col), !!sym(target_col))
    
    return(biodict)
}



#' Convert genes between species efficiently
#'
#' @param genes Vector of gene symbols to convert
#' @param org.from Source organism ("human" or "mouse")
#' @param org.to Target organism ("human" or "mouse")
#' @param one.to.many If TRUE, returns all possible mappings
#' @return Vector of converted gene symbols
#' @importFrom dplyr %>% filter select group_by slice distinct inner_join
#' @export
convert_genes <- function(genes, org.from = "human", org.to = "mouse", one.to.many = FALSE) {
    # Validate inputs
    validate_org(org.from, genes = genes)
    validate_org(org.to)
    stopifnot(org.from != org.to)
    
    # Get conversion dictionary
    biodict <- get_conversion_dict(org.from, org.to)
    source_col <- paste0(org.from, "_gene_symbol")
    target_col <- paste0(org.to, "_gene_symbol")
    
    # Convert input genes to data frame
    genes_df <- data.frame(source_gene = genes)
    
    # Join with dictionary
    result <- genes_df %>%
        dplyr::inner_join(biodict, by = c("source_gene" = source_col))
    
    # Process results based on one.to.many parameter
    if (one.to.many) {
        # Return all unique target genes
        unique(result[[target_col]])
    } else {
        # Get unique source genes and their first target
        result %>%
            dplyr::group_by(source_gene) %>%
            dplyr::slice(1) %>%
            dplyr::pull(!!sym(target_col))
    }
}



#' Convert data frame between species efficiently
#'
#' @param df Data frame containing gene symbols
#' @param gene_column Column name containing gene symbols
#' @param org.from Source organism ("human" or "mouse")
#' @param org.to Target organism ("human" or "mouse")
#' @param one.to.many If TRUE, returns all possible mappings
#' @return Data frame with converted gene symbols
#' @importFrom dplyr %>% filter group_by summarize select mutate left_join inner_join distinct
#' @export
convert_df <- function(df, gene_column, org.from = "human", org.to = "mouse", one.to.many = FALSE) {
    # Validate inputs
    df <- validate_df(df, gene_column)
    validate_org(org.from, genes = df[[gene_column]])
    validate_org(org.to)
    stopifnot(org.from != org.to)
    
    # Get conversion dictionary
    biodict <- get_conversion_dict(org.from, org.to)
    source_col <- paste0(org.from, "_gene_symbol")
    target_col <- paste0(org.to, "_gene_symbol")
    
    # Apply grouping and summarization
    dict_summary <- biodict %>%
        dplyr::group_by(!!sym(source_col)) %>%
        dplyr::summarize(
            !!sym(target_col) := paste(!!sym(target_col), collapse = ", ")
        )
    
    # Merge with original data frame
    converted_df <- df %>%
        dplyr::left_join(dict_summary, 
                         by = setNames(source_col, gene_column))
    
    # Process one-to-many relationships if needed
    if (!one.to.many) {
        converted_df[[target_col]] <- gsub(",.*", "", converted_df[[target_col]])
    }
    
    return(converted_df)
}
   
#' Convert expression matrix or data frame between species efficiently
#'
#' @param exprs Expression matrix, data frame, or data.table
#' @param org.from Source organism ("human" or "mouse")
#' @param org.to Target organism ("human" or "mouse")
#' @param many.to.one If TRUE, aggregates multiple mappings
#' @param normalized If TRUE, uses mean; if FALSE, uses sum
#' @return Converted expression data with same type as input
#' @importFrom dplyr %>% filter group_by summarize select mutate left_join inner_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
convert_exprs <- function(exprs, org.from = "human", org.to = "mouse", 
                         many.to.one = TRUE, normalized = FALSE) {
    
    # Validate inputs
    exprs <- validate_exprs(exprs)
    validate_org(org.from, genes = rownames(exprs))
    validate_org(org.to)
    stopifnot(org.from != org.to)
    
    # Get conversion dictionary and column names
    biodict <- get_conversion_dict(org.from, org.to)
    source_col <- paste0(org.from, "_gene_symbol")
    target_col <- paste0(org.to, "_gene_symbol")
    
    if(many.to.one){
        converted_exprs <- exprs %>%
            rownames_to_column("gene") %>%
            pivot_longer(!gene, names_to = "sample", values_to = "exprs") %>%
            merge(., biodict %>% distinct(!!sym(source_col), .keep_all = TRUE), by.x = "gene", by.y = source_col, all.x = TRUE) %>%
            filter(!!sym(target_col) != "NA") %>%
            pivot_wider(names_from = sample, values_from = exprs) %>%
            dplyr::select(-c(gene))

        converted_vector <- converted_exprs[[target_col]]
        converted_exprs[[target_col]] <- NULL
        converted_exprs <- summarize_genes(converted_exprs, gene_sym_vec = converted_vector, normalized = normalized)}
    else{
        dict_summary <- biodict %>% 
            group_by(!!sym(source_col)) %>%
            mutate(hn = n()) %>%
            group_by(!!sym(target_col)) %>%
            mutate(mn = n()) %>%
            filter(hn == 1 & mn == 1)
        dict_summary <- dict_summary[c(1,2)]
        converted_exprs <- exprs[which(rownames(exprs) %in% dict_summary[[source_col]]),]

        converted_exprs <- converted_exprs %>%
            merge(., dict_summary, by.x = 0, by.y = source_col, all.x = TRUE) %>%
            dplyr::select(-c("Row.names")) %>%
            column_to_rownames(target_col)}
    return(converted_exprs)
}