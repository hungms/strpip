#' Summarize Gene Expression
#'
#' This function aggregates expression values from multiple isoforms of the same gene
#' into a single gene-level expression value.
#'
#' @param input Gene expression data with isoforms/peaks in rows and samples in columns.
#'              Can be a data frame, data.table, or matrix.
#' @param gene_sym_vec Vector of gene symbols corresponding to each isoform/peak
#' @param normalized Logical. If TRUE, calculates the mean expression across isoforms/peaks.
#'                   If FALSE, calculates the sum of expression across isoforms/peaks.
#'                   Default is FALSE
#' @return A data frame, data.table, or matrix with unique gene names as row names.
#'         Output type matches the class of the input.
#' @importFrom data.table as.data.table setDT copy setkey melt dcast :=
#' @export
summarize_genes <- function(input, gene_sym_vec, normalized = FALSE) {
    # Validate input
    stopifnot(
        "Input cannot be empty" = nrow(input) > 0 && ncol(input) > 0,
        "gene_sym_vec must have same length as number of rows" = length(gene_sym_vec) == nrow(input),
        "gene_sym_vec cannot contain NA values" = !any(is.na(gene_sym_vec))
    )
    
    # Store input properties
    is_matrix <- is.matrix(input)
    is_dt <- is.data.table(input)
    orig_colnames <- colnames(input)
    
    # Convert to data.table efficiently
    dt <- if (is_matrix) {
        setDT(as.data.table(input))
    } else if (is_dt) {
        copy(input)
    } else {
        setDT(as.data.table(input))
    }
    
    # Add gene names efficiently
    dt[, Gene.Name := gene_sym_vec]
    
    # Set key for faster aggregation
    setkey(dt, Gene.Name)
    
    # Aggregate efficiently
    if (normalized) {
        result <- dt[, lapply(.SD, mean, na.rm = TRUE), by = Gene.Name]
    } else {
        result <- dt[, lapply(.SD, sum, na.rm = TRUE), by = Gene.Name]
    }
    
    # Process result efficiently
    gene_col <- result$Gene.Name
    result[, Gene.Name := NULL]
    
    # Return in original format
    if (is_matrix) {
        result_mat <- as.matrix(result)
        if (!is.null(orig_colnames)) {
            result_mat <- result_mat[, orig_colnames]
        }
        rownames(result_mat) <- gene_col
        result_mat
    } else if (is_dt) {
        setattr(result, "row.names", gene_col)
        result[]
    } else {
        setattr(as.data.frame(result), "row.names", gene_col)
    }
}