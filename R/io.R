#' Convert List to Data Frame
#'
#' Converts a list of vectors into a data frame, padding shorter vectors with empty strings
#' to match the length of the longest vector.
#'
#' @param list A list of vectors to convert
#' @return A data frame where each column corresponds to a vector from the input list.
#'         Shorter vectors are padded with empty strings to match the length of the longest vector.
#' @examples
#' \dontrun{
#' my_list <- list(c("A", "B", "C"), c("D", "E"), c("F", "G", "H", "I"))
#' result <- convert_list_to_df(my_list)
#' print(result)
#' }
#' @export
convert_list_to_df <- function(list) {
    # Validate input
    if (!is.list(list)) {
        stop("Input must be a list", call. = FALSE)
    }
    
    # Check if list is empty
    if (length(list) == 0) {
        return(data.frame())
    }
    
    # Pad shorter vectors with empty strings
    max_length <- max(sapply(list, length))
    padded.list <- lapply(list, function(v) {
        c(v, rep("", max_length - length(v)))
    })
    
    # Convert to data frame
    df <- as.data.frame(padded.list, stringsAsFactors = FALSE)
    return(df)
}

#' Read GMT File
#'
#' Reads a Gene Matrix Transposed (GMT) file and converts it to a data frame.
#' GMT files are commonly used for gene set collections.
#'
#' @param gmt Path to the GMT file
#' @return A data frame where:
#'         - Rows are genes
#'         - Columns are gene sets
#'         - Values indicate gene set membership
#' @examples
#' \dontrun{
#' gmt_data <- read_gmt("path/to/genesets.gmt")
#' head(gmt_data)
#' }
#' @export
read_gmt <- function(gmt) {
    # Check if file exists
    if (!file.exists(gmt)) {
        stop("GMT file not found: ", gmt, call. = FALSE)
    }
    
    # Try to read and process the file
    tryCatch({
        gmt.df <- utils::read.table(gmt, sep = "\t", header = FALSE, row.names = 1)
        
        # Check if file has expected format
        if (ncol(gmt.df) < 1) {
            stop("Invalid GMT file format: file must have at least two columns", call. = FALSE)
        }
        
        # Process the data
        gmt.df <- gmt.df[, -1] %>% 
                  t() %>% 
                  as.data.frame()
                  
        return(gmt.df)
    }, error = function(e) {
        # Handle read errors
        stop("Error reading GMT file: ", e$message, call. = FALSE)
    })
}

#' Read GCT File
#'
#' Reads a Gene Cluster Text (GCT) file and converts it to a data frame.
#' GCT files are commonly used for gene expression data.
#'
#' @param gct Path to the GCT file
#' @return A data frame where:
#'         - Rows are genes
#'         - Columns are samples
#'         - Values are expression measurements
#' @examples
#' \dontrun{
#' gct_data <- read_gct("path/to/expression.gct")
#' head(gct_data)
#' }
#' @export
read_gct <- function(gct) {
    # Check if file exists
    if (!file.exists(gct)) {
        stop("GCT file not found: ", gct, call. = FALSE)
    }
    
    # Try to read and process the file
    tryCatch({
        # Read GCT file
        gct.df <- utils::read.table(gct, sep = "\t", header = TRUE, row.names = 1, skip = 2)
        
        # Check if file has expected format
        if (ncol(gct.df) < 1) {
            stop("Invalid GCT file format", call. = FALSE)
        }
        
        # Remove Description column if present
        if ("Description" %in% colnames(gct.df)) {
            gct.df <- gct.df[, -which(colnames(gct.df) == "Description")]
        } else {
            gct.df <- gct.df[, -1]
        }
        
        return(gct.df)
    }, error = function(e) {
        # Handle read errors
        stop("Error reading GCT file: ", e$message, call. = FALSE)
    })
}

#' Write GMT File
#'
#' Saves a data frame as a Gene Matrix Transposed (GMT) file.
#' The output file will be formatted according to the GMT specification.
#'
#' @param df Data frame to save
#' @param save_name Name of the output file
#' @param save_dir Directory to save the file in. Defaults to current working directory.
#' @return No return value. Creates a GMT file at the specified location.
#' @examples
#' \dontrun{
#' write_gmt(gene_sets_df, "my_genesets.gmt", "output/")
#' }
#' @export
write_gmt <- function(df, save_name = "", save_dir = getwd()) {
    # Validate inputs
    if (!is.data.frame(df)) {
        stop("Input must be a data frame", call. = FALSE)
    }
    
    if (save_name == "") {
        stop("Output filename cannot be empty", call. = FALSE)
    }
    
    if (!dir.exists(save_dir)) {
        stop("Output directory does not exist: ", save_dir, call. = FALSE)
    }
    
    # Create output path
    output_path <- file.path(save_dir, save_name)
    
    # Transform and write data
    tryCatch({
        # Transform data frame to GMT format
        gmt.df <- df %>% 
            t() %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column("Name")
            
        # Add description column (usually same as gene set name)
        gmt.df$Description <- gmt.df$Name
        
        # Reorder columns
        gmt.df <- gmt.df[, c(ncol(gmt.df), 1:(ncol(gmt.df)-1))]
        
        # Write to file
        utils::write.table(gmt.df, file = output_path, 
                          col.names = FALSE, row.names = FALSE, 
                          quote = FALSE, sep = "\t")
                          
        # Success message
        invisible(output_path)
    }, error = function(e) {
        stop("Error writing GMT file: ", e$message, call. = FALSE)
    })
}

#' Write GCT File
#'
#' Saves a data frame as a Gene Cluster Text (GCT) file.
#' The output file will be formatted according to the GCT specification.
#'
#' @param df Data frame to save
#' @param save_name Name of the output file
#' @param save_dir Directory to save the file in. Defaults to current working directory.
#' @return No return value. Creates a GCT file at the specified location.
#' @examples
#' \dontrun{
#' write_gct(expression_df, "my_expression.gct", "output/")
#' }
#' @export
write_gct <- function(df, save_name = "", save_dir = getwd()) {
    # Validate inputs
    if (!is.data.frame(df)) {
        stop("Input must be a data frame", call. = FALSE)
    }
    
    if (save_name == "") {
        stop("Output filename cannot be empty", call. = FALSE)
    }
    
    if (!dir.exists(save_dir)) {
        stop("Output directory does not exist: ", save_dir, call. = FALSE)
    }
    
    # Create output path
    output_path <- file.path(save_dir, save_name)
    
    # Transform and write data
    tryCatch({
        # Create GCT format matrix
        gct.df <- matrix(ncol = ncol(df) + 2, nrow = nrow(df) + 3)
        
        # Add header information
        gct.df[1, 1] <- "#1.3"
        gct.df[2, 1:2] <- dim(df)
        gct.df[3, 1:2] <- c("Name", "Description")
        gct.df[3, 3:ncol(gct.df)] <- colnames(df)
        
        # Add row names and data
        gct.df[4:(nrow(gct.df)), 1] <- rownames(df)
        gct.df[4:(nrow(gct.df)), 2] <- rownames(df)  # Description = Name
        gct.df[4:nrow(gct.df), 3:ncol(gct.df)] <- as.matrix(df)
        
        # Convert to data frame
        gct.df <- as.data.frame(gct.df)
        
        # Write to file
        utils::write.table(gct.df, file = output_path, 
                          col.names = FALSE, row.names = FALSE, 
                          quote = FALSE, sep = "\t")
                          
        # Success message
        invisible(output_path)
    }, error = function(e) {
        stop("Error writing GCT file: ", e$message, call. = FALSE)
    })
}

#' Summarize Gene Expression
#'
#' This function aggregates expression values from multiple isoforms of the same gene
#' into a single gene-level expression value. It can either sum or average the isoform
#' expression values for each gene.
#'
#' @param df Gene expression dataframe with isoforms/peaks in rows and samples in columns
#' @param gene_sym_vec Vector of gene symbols corresponding to each isoform/peak
#' @param normalized Logical. If TRUE, calculates the mean expression across isoforms/peaks.
#'                   If FALSE, calculates the sum of expression across isoforms/peaks.
#'                   Default is FALSE
#' @return A gene expression dataframe with unique gene names as row names
#' @examples
#' # Create example data
#' isoform_expr <- data.frame(
#'   sample1 = c(10, 20, 15, 30),
#'   sample2 = c(12, 25, 18, 35),
#'   row.names = c("isoform1", "isoform2", "isoform3", "isoform4")
#' )
#' gene_symbols <- c("geneA", "geneA", "geneB", "geneB")
#' 
#' # Sum isoform/peak expression
#' gene_expr <- summarize_genes(isoform_expr, gene_symbols)
#' 
#' # Average isoform/peak expression
#' gene_expr <- summarize_genes(isoform_expr, gene_symbols, normalized = TRUE)
#' @export
summarize_genes <- function(df, gene_sym_vec, normalized = FALSE) {
    # Validate inputs
    if (!is.data.frame(df)) {
        stop("Expression data must be a data frame", call. = FALSE)
    }
    
    if (!is.character(gene_sym_vec)) {
        stop("Gene symbols must be a character vector", call. = FALSE)
    }
    
    if (nrow(df) != length(gene_sym_vec)) {
        stop("Number of rows in data frame must match length of gene symbol vector", call. = FALSE)
    }
    
    # Transform and summarize data
    tryCatch({
        # Add gene names to data frame
        df$Gene.Name <- gene_sym_vec
        
        # Reshape to long format for easier aggregation
        df_long <- df %>%
            as.data.frame() %>%
            tidyr::pivot_longer(!c("Gene.Name"), names_to = "samples", values_to = "exprs") %>%
            dplyr::group_by(.data$Gene.Name, .data$samples)
    
        # Calculate summary statistic based on normalized parameter
        if (normalized) {
            df_summary <- df_long %>% dplyr::summarize(exprs = mean(.data$exprs), .groups = "drop")
        } else {
            df_summary <- df_long %>% dplyr::summarize(exprs = sum(.data$exprs), .groups = "drop")
        }
    
        # Convert back to wide format
        df_result <- df_summary %>% 
            tidyr::pivot_wider(names_from = "samples", values_from = "exprs") %>%
            tibble::column_to_rownames("Gene.Name")
            
        return(df_result)
    }, error = function(e) {
        stop("Error summarizing gene expression: ", e$message, call. = FALSE)
    })
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL 