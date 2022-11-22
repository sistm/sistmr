#' Color the dendrogram labels
#' 
#' Assign colors to the labels of a \code{hclust} dendrogram according to a factor 
#'
#' @param hclust a \code{hclust} object.
#' @param labels a character vector of labels for the leaves of the tree. By default, the labels in the \code{hclust} object are used if specified.
#' @param lab.col a color vector for the labels. 
#' @param ... further arguments passed to plot function
#'
#' @return a dendrogram
#' 
#' @importFrom dendextend labels_colors<- is.hclust
#' @importFrom graphics plot.new
#' @importFrom grDevices recordPlot
#' @importFrom stats as.dendrogram order.dendrogram
#' 
#' @author Mélanie Huchon
#' 
#' @export
#'
#' @examples
#' #Load data example
#' data(iris)
#' #Run hclust on this data
#' hc <- hclust(dist(iris[, -5]))
#' #Color vector 
#' color <- c("#440154FF", "#21908CFF", "#FDE725FF")
#' #Use function to color labels
#' colorDendro(hc, labels = iris[,5], lab.col = color, main = "Dendrogram of Iris data")

colorDendro <- function(hclust, labels = hclust$labels, lab.col = NULL, ...){
  
  #CHECKS :
  #on class of hclust
  if(!is.hclust(hclust)){
    stop("*hclust* sould be a *hclust* object")
  }
  #on labels
  if(is.null(labels)) {
    stop("Labels don't found in *hclust* object, please specify a character vector of labels.")
    q("no", status = 1, runLast = FALSE)
  }
  if(is.null(hclust$labels) | !identical(hclust$labels, labels)){
    hclust$labels <- labels
  }
  #on length of color and label vectors
  if(!is.null(lab.col) & !identical(length(unique(labels)), length(lab.col))){
    stop("Number of labels and colors don't correspond.")
  }
  if(is.null(lab.col)){
    lab.col <- 1:length(unique(labels))
  }
  
  #Association between labels and colors
  tab_labelCol <- data.frame(labels = factor(labels))
  for(i in 1:length(unique(labels))){
    tab_labelCol[tab_labelCol$labels == levels(tab_labelCol$labels)[i], "color"] <- unique(lab.col)[i]
  }
  
  #Calculate the dendrogram
  hc <- as.dendrogram(hclust)
  
  #Retrieve colors and sort them based on their order in hc
  colors_to_use <- as.character(tab_labelCol$color[order.dendrogram(hc)])

  #Apply color on dendrogram
  labels_colors(hc) <- colors_to_use

  # Plotting dendrogram
  plot(hc, ...)
}
