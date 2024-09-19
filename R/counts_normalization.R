#' Calulate CPM Normalization for Raw Gene Expression Counts 
#'
#' @param raw_counts numeric matrix or data frame of raw gene expression counts.  
#' @param MARGIN "row" or "col" to identify the place of sample (library size). If the \code{raw_counts} matrix has a \code{nxG} size, put "row" if it's opposite put "col".
#' @param log2 a logical indicating whether the normalized counts will be transformed in log2. Default is \code{TRUE}.
#' @param plot a logical indicating whether a graph will be displayed for the normalization check. Default is \code{TRUE}.
#' @param TMM_method a logical indicating whether the library size will be normalized by TMM method of \code{edgeR} package. Default is \code{FALSE}.
#'
#' @return normalized counts matrix 
#' 
#' @importFrom utils menu
#' @importFrom edgeR DGEList calcNormFactors
#'  
#' @author Mélanie Huchon
#' 
#' @export
#' 
#' @examples
#' 
#' #Generate data
#' sim_Genes <- c(runif(90,0,1), rbinom(60,2,0.4), runif(30,1250,2000), runif(20,100,300))
#' raw_data <- matrix(sim_Genes, nrow = 10)
#' colnames(raw_data) <- paste0("G", sample(1:ncol(raw_data), ncol(raw_data), replace = FALSE))
#' rownames(raw_data) <- paste0("Sample_", 1:nrow(raw_data))
#' 
#' #Use the method
#' norm_data <- counts_normalization(raw_counts = data.frame(raw_data), MARGIN = "row")

counts_normalization <- function(raw_counts, MARGIN = c("row", "col"), log2 = TRUE, plot = TRUE, TMM_method = FALSE){
  
  # Ensure that MARGIN is either "row" or "col"
  MARGIN <- match.arg(MARGIN) 
  
  if (MARGIN == "row" && ncol(raw_counts) < nrow(raw_counts)) {
    input <- menu(c("Yes", "No"), title = "Are you sure that your samples are in rows? (nb_col < nb_row)")
    if (input == 2) stop("The normalization has been stopped. Please check your input dimensions.")
  }
  
  if (MARGIN == "col" && ncol(raw_counts) > nrow(raw_counts)) {
    input <- menu(c("Yes", "No"), title = "Are you sure that your samples are in columns? (nb_col > nb_row)")
    if (input == 2) stop("The normalization has been stopped. Please check your input dimensions.")
  }
  
  # Transpose if samples are in rows
  if (MARGIN == "row") {
    raw_counts <- t(raw_counts)
  }
  
  # Calculate library sizes (sum of counts for each sample)
  lib.size <- colSums(raw_counts)
  
  # Apply TMM normalization if requested
  if (TMM_method) {
    # Create DGEList object and calculate TMM normalization factors
    dge <- DGEList(counts = raw_counts)
    dge <- calcNormFactors(dge, method = "TMM")
    # Adjust library sizes with TMM normalization factors
    lib.size <- lib.size * dge$samples$norm.factors
  }
  # Normalize counts by library size (counts per million)
  norm_counts <- sweep(raw_counts, 2, lib.size, "/") * 1e6
  
  # Transpose back if samples were originally in rows
  if (MARGIN == "row") {
    norm_counts <- t(norm_counts)
  }
  
  # Apply log2 transformation if requested
  if (log2) {
    norm_counts <- log2(norm_counts + 1)
  }
  
  # Plot normalized counts if requested
  if (plot) {
    if (log2) {
      norm_counts_plot <- 2^norm_counts
    } else {
      norm_counts_plot <- norm_counts
    }
    
    barplot(if (MARGIN == "row") rowSums(norm_counts) else colSums(norm_counts), 
            xlab = "Samples", main = "Normalized Counts")
  }
  
  return(as.data.frame(norm_counts))
}
