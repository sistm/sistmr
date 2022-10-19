#' Multiple boxplots for many times
#' 
#' @param data a dataset from which the variable \code{x_var} and \code{y_var} should be taken.
#' @param x_var corresponding to the x coordinates for the plot, it must be a factor to obtain multiple boxplots.
#' @param y_var corresponding to the y coordinates for the plot.
#' @param add_points if you want to add points on boxplots. Default value is \code{TRUE}. 
#'
#' @return a \code{ggplot2} object
#' @import ggbeeswarm
#' @import ggplot2
#' @import rlang 
#' @export
#'
#' @examples 
#' library(ggplot2)
#' 
#' #Generate data
#' x_ex <- factor(c(rep("J0", 10), rep("J7", 10), rep("J14", 10)), levels = c("J0", "J7", "J14"))
#' y_ex <- rnorm(30)
#' 
#' data_ex <- cbind.data.frame(x_ex, y_ex)
#' 
#' #Plotting
#' multipleBoxplots(data = data_ex, x_var = x_ex, y_var = y_ex)
#' 
#' multipleBoxplots(data = data_ex, x_var = x_ex, y_var = y_ex) + 
#' labs(x = "Time", y = "Value") + 
#' theme(legend.position = "none")

multipleBoxplots <- function(data, x_var, y_var, add_points = TRUE){

  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  x_var <- enquo(x_var)
  y_var <- enquo(y_var)
  
  x_val <- eval_tidy(x_var, data)
  if(!is.factor(x_val)){
    stop("x_var must be a factor!")
  }
  
  #numbers of factors
  nb_factors <- length(levels(x_val))
  
  #Initial plot 
  plot <- ggplot(data) +
          #To add boxplots
          geom_boxplot(aes(y = !!y_var, x = !!x_var, color = !!x_var, fill = !!x_var), outlier.shape = NA, width = 0.4, lwd = 0.6) + 
          #To see points and don't have background color in boxplots
          scale_fill_manual(values = rep("transparent", nb_factors)) + 
          #To don't have background color in legend
          guides(fill = "none") +
          #To change background plot
          theme_classic()
  
  if(add_points){
    #To add points on graph
    plot <- plot + geom_quasirandom(aes(y = !!y_var, x = !!x_var, color = !!x_var), size = 1, alpha = 0.5) 
  }

  return(plot)
}

