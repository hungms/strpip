
#' remove_xy_genes
#'
#' Remove gene symbols from the XY chromosome from dataframe.
#' @param df Dataframe.
#' @param gene_col Column name of gene symbols.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @return Dataframe with no XY genes.
#' @export
remove_xy_genes <- function(df, gene_col = "rownames(.)", org = "human"){
    df <- df %>% dplyr::filter(!!sym(gene_col) %in% get_xy_genes(org = org))
    return(df)}

#' remove_mt_genes
#'
#' Remove gene symbols from the MT chromosome from dataframe.
#' @param df Dataframe.
#' @param gene_col Column name of gene symbols.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @return Dataframe with no MT genes.
#' @export
remove_mt_genes <- function(df, gene_col = "rownames(.)", org = "human"){
    df <- df %>% dplyr::filter(!!sym(gene_col) %in% get_mt_genes(org = org))
    return(df)}

#' remove_vdj_genes
#'
#' Remove gene symbols from the MT chromosome from dataframe.
#' @param df Dataframe.
#' @param gene_col Column name of gene symbols.
#' @param org Organism either "human" or "mouse". Default is "human".
#' @return Dataframe with no MT genes.
#' @export
remove_vdj_genes <- function(df, gene_col = "rownames(.)"){
    df <- df %>% dplyr::filter(!!sym(gene_col) %in% )
    return(df)}

#' remove_genes_by_string
#'
#' Remove gene symbols matching a string from dataframe.
#' @param df Dataframe.
#' @param gene_col Column name of gene symbols.
#' @return Dataframe with genes removed.
#' @export
remove_genes <- function(df, gene_col = "rownames(.)", string = ""){
    df <- df %>% dplyr::filter(!!sym(gene_col) %in% get_mt_genes(org = org))
    return(df)}




