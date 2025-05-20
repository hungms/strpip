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
#' result <- list_to_df(my_list)
#' print(result)
#' }
#' @export
list_to_df <- function(list) {
    # Validate input
    stopifnot("Input must be a list" = is.list(list))
    
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

#' Convert Data Frame to List
#'
#' Converts a data frame into a list of vectors, removing empty strings.
#'
#' @param df A data frame to convert
#' @return A list of vectors where each element corresponds to a column from the input data frame.
#' @examples
#' \dontrun{
#' my_df <- data.frame(c("A", "B", "C"), c("D", "E"), c("F", "G", "H", "I"))
#' result <- df_to_list(my_df)
#' print(result)
#' }
#' @export
df_to_list <- function(df) {
    # Validate input
    stopifnot("Input must be a data frame" = is.data.frame(df))
    list <- lapply(df, function(x) x[!is.na(x)])
    return(list)}

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
  stopifnot(file.exists(gmt))
  lines <- readLines(gmt)
  split_lines <- strsplit(lines, "\t")
  max_cols <- max(sapply(split_lines, length))
  padded <- lapply(split_lines, function(x) { length(x) <- max_cols; x })
  df <- data.frame(do.call(rbind, padded), stringsAsFactors = FALSE) %>% t()
  colnames(df) <- df[1,]
  df <- df[-c(1,2),]
  rownames(df) <- NULL
  df <- as.data.frame(df)
  return(df)
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
    stopifnot("GCT file not found" = file.exists(gct))
    
    # Read GCT file
    gct.df <- utils::read.table(gct, sep = "\t", header = TRUE, row.names = 1, skip = 2)
    
    # Check if file has expected format
    stopifnot("Invalid GCT file format" = ncol(gct.df) >= 1)
    
    # Remove Description column if present
    if ("Description" %in% colnames(gct.df)) {
        gct.df <- gct.df[, -which(colnames(gct.df) == "Description")]} 
    else {
        gct.df <- gct.df[, -1]}
    
    return(gct.df)
}

#' Write GMT File
#'
#' Saves a data frame or list as a Gene Matrix Transposed (GMT) file.
#' The output file will be formatted according to the GMT specification.
#'
#' @param input A data frame where columns are gene sets and values are genes,
#'              or a named list where each element is a character vector of genes
#' @param file Output file path
#' @return No return value. Creates a GMT file at the specified location.
#' @examples
#' \dontrun{
#' # Using a data frame
#' write_gmt(gene_sets_df, "my_genesets.gmt")
#' 
#' # Using a list
#' gene_sets_list <- list(
#'   pathway1 = c("gene1", "gene2", "gene3"),
#'   pathway2 = c("gene2", "gene4", "gene5")
#' )
#' write_gmt(gene_sets_list, "my_genesets.gmt")
#' }
#' @export
write_gmt <- function(input, file) {
  
  # Validate file path
  stopifnot("File path must be a character string" = is.character(file))
  stopifnot("Directory does not exist" = file.exists(dirname(file)))
  
  # Convert list to data frame if needed
  if (is.list(input) && !is.data.frame(input)) {
    stopifnot("List must have names" = !is.null(names(input)))
    stopifnot("List must contain character vectors" = all(sapply(input, is.character)))
    
    # Convert list to data frame using list_to_df
    input <- list_to_df(input)
  }
  
  # Validate that input is now a data frame
  stopifnot("Input must be a data frame or a named list of character vectors" = is.data.frame(input))
  
  # Open file connection
  con <- file(file, open = "wt")
  
  # Process each gene set (column)
  for (gene_set_name in colnames(input)) {
    genes <- na.omit(input[[gene_set_name]])
    genes <- genes[genes != ""]
    
    # Prepare the line: gene set name, description (same as name), genes
    line <- c(gene_set_name, gene_set_name, genes)
    writeLines(paste(line, collapse = "\t"), con)
  }
  
  # Close connection
  close(con)
  cat("GMT file written to", file, "\n")
  
  invisible(file)
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
write_gct <- function(df, file) {
    # Validate inputs
    stopifnot(
        "Input must be a data frame" = is.data.frame(df),
        "Output filename cannot be empty" = basename(file) != "",
        "Output directory does not exist" = dir.exists(dirname(file)))
    
    # Create output path
    output_path <- file
    
    # Create GCT format matrix
    gct.df <- matrix("", ncol = ncol(df) + 2, nrow = nrow(df) + 3)
    
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
    utils::write.table(gct.df, file = file, 
                      col.names = FALSE, row.names = FALSE, 
                      quote = FALSE, sep = "\t")
                      
    # Success message
    invisible(file)
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