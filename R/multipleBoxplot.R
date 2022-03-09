#' Multiple boxplots for many times
#' 
#' @param data a dataset from which the variable \code{x_var} and \code{y_var} should be taken.
#' @param x_var corresponding to the x coordinates for the plot.
#' @param y_var corresponding to the y coordinates for the plot.
#' 
#' @return a \code{ggplot2} object
#' @import ggbeeswarm 
#' @import ggplot2
#' @export
#'

multipleBoxplots <- function(data, x_var, y_var){

  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  x_var <- enquo(x_var)
  y_var <- enquo(y_var)
  
  # browser()
  
  #Initial plot 
  plot <- ggplot(data) +
          #To add points on graph
          geom_quasirandom(aes(y = !!y_var, x = !!x_var, color = !!x_var), size = 1, alpha = 0.5) +
          #To add boxplots
          geom_boxplot(aes(y = !!y_var, x = !!x_var, color = !!x_var, fill = !!x_var), outlier.shape = NA, width = 0.4, lwd = 0.6) + 
          #To see points and don't have background color in boxplots
          scale_fill_manual(values = c("transparent", "transparent", "transparent")) +
          #To don't have background color in legend
          guides(fill = "none") +
          #To change background plot
          theme_classic()

  return(plot)
}

