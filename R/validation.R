#' Validate Org Parameters
#'
#' @param org Organism type, either "human" or "mouse"
#' @param genes Gene symbols to validate
#' @return Logical indicating whether the genes are valid
#' @keywords internal
validate_org <- function(org, genes = NULL) {
    stopifnot("org must be either 'human' or 'mouse'" = org %in% c("human", "mouse"))
    if (!is.null(genes)) {
        if (org == "mouse") {
            stopifnot(any(grepl("[a-z]", genes))) } 
        if (org == "human") {
            stopifnot(!all(grepl("[a-z]", genes))) } }
}

#' Validate input data structure
#'
#' @param df A data frame or matrix to validate
#' @param gene_column Column name containing gene symbols (for data frames) or NULL for matrix
#' @return Validated data frame
#' @importFrom dplyr %>% filter mutate select
#' @keywords internal
validate_df <- function(df, gene_column = NULL) {

    stopifnot(is.data.frame(df) || is.matrix(df))

    # convert matrix to data frame
    if (is.matrix(df)) {
        genes <- rownames(df)
        df <- as.data.frame(df)
        rownames(df) <- genes}
    
    # check for gene column
    if (!is.null(gene_column)) {
        if (!gene_column %in% names(df)) {
            stop("Gene column '", gene_column, "' not found in data frame")}
        
        # filter for valid gene symbols
        orig_rows <- nrow(df)
        df <- df %>%
            dplyr::filter(!is.na(!!sym(gene_column))) %>%
            dplyr::filter(grepl("[A-Z]", !!sym(gene_column)))
        
        # Report changes
        removed <- orig_rows - nrow(df)
        if (removed > 0) {
            message("Removed ", removed, " rows with NA values or invalid gene symbols")}
    }
    
    return(df)
}

#' Validate expression data
#' 
#' @param exprs Expression matrix or data frame to validate
#' @return Validated data frame with gene column
#' @importFrom dplyr %>% filter select mutate
#' @keywords internal
validate_exprs <- function(exprs) {

    # check input type
    stopifnot(is.data.frame(exprs) || is.matrix(exprs))

    # Convert to data frame if needed
    if (is.matrix(exprs)) {
        genes <- rownames(exprs)
        exprs <- as.data.frame(exprs)
        rownames(exprs) <- genes}
    
    # Filter invalid rows
    valid_idx <- !is.na(rownames(exprs)) & 
                 !grepl("^NA(\\.[0-9]+)?$", rownames(exprs)) & 
                 grepl("[A-Z]", rownames(exprs))
    
    # Report changes
    removed <- sum(!valid_idx)
    if (removed > 0) {
        message("Removed ", removed, " rows with NA values or invalid gene symbols") }
    
    exprs <- exprs[valid_idx, ]
    return(exprs)
}