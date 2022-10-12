#' Calulate CPM Normalization for Raw Gene Expression Counts 
#'
#' @param raw_counts numeric matrix or data frame of raw gene expression counts.  
#' @param MARGIN "row" or "col" to identify the place of sample (library size). If the \code{raw_counts} matrix has a \code{nxG} size, put "row" if it's opposite put "col".
#' @param log2 a logical indicating whether the normalized counts will be transformed in log2. Default is \code{TRUE}.
#' @param plot a logical indicating whether a graph will be displayed for the normalization check. Default is \code{TRUE}.
#'
#' @return normalized counts matrix 
#' 
#' @importFrom utils menu
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

counts_normalization <- function(raw_counts, MARGIN = c("row", "col"), log2 = TRUE, plot = TRUE){
  
  input = NULL

  # Check if the dimensions and MARGIN condition correspond and ask if the normalisation is on the dimension 
  if(ncol(raw_counts) < nrow(raw_counts) & identical(MARGIN, "row")){
    input <- menu(c("Yes", "No"), title = "Are you sure that your samples are in rows? (nb_col < nb_row)") 
    if(input == 1){
      norm_counts <- apply(raw_counts, MARGIN = 1, function(v) {
        (v + 0.5)/(sum(v) + 1) * 10^6
      })
    }
  }
  
  if(ncol(raw_counts) > nrow(raw_counts) & identical(MARGIN, "col")){
    input <- menu(c("Yes", "No"), title = "Are you sure that your samples are in columns? (nb_col > nb_row)") 
    if(input == 1){
      norm_counts <- apply(raw_counts, MARGIN = 2, function(v) {
        (v + 0.5)/(sum(v) + 1) * 10^6
      })
    }
  }
  
  #If the answer is "No", the script stop.
  if(identical(input, 2)){
    opt <- options(show.error.messages=FALSE, browser = FALSE)
    on.exit(options(opt))
  #If the conditions are ok, calculate the normalized counts.
  }else{
    if(identical(MARGIN, "row")){
      norm_counts <- apply(raw_counts, MARGIN = 1, function(v) {
        (v + 0.5)/(sum(v) + 1) * 10^6
      })
    }else if(identical(MARGIN, "col")){
      norm_counts <- apply(raw_counts, MARGIN = 2, function(v) {
        (v + 0.5)/(sum(v) + 1) * 10^6
      })
    }
  }
  #To transform in log2 
  if(log2){
    norm_counts <- log2(norm_counts)
  }
  #To transpose if necessary
  if(!identical(dim(raw_counts)[1], dim(norm_counts)[1])){
    norm_counts <- t(norm_counts)
  }
  
  # To display the graph control
  if(plot){
    if(log2){
      norm_counts_plot <- 2^norm_counts
    }else{
      norm_counts_plot <- norm_counts
    }
    if(identical(MARGIN, "row")){
      barplot(rowSums(norm_counts_plot), xlab = "Samples")
    }else if(identical(MARGIN, "col")){
      barplot(colSums(norm_counts_plot), xlab = "Samples")
    }
  }

  return(as.data.frame(norm_counts))
}
