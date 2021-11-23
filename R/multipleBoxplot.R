#' Multiple boxplots for many times
#' 
#' @param abs corresponding to the x coordinates for the plot.
#' @param ord corresponding to the y coordinates for the plot.
#' @param data Dataset to use for plot.
#' @param color_palette a color palette for boxplot.
#' @param color_manual a color vector where each color is for one boxplot.
#'
#' @return 
#' @import ggplot2
#' @import ggbeeswarm
#' @export

multipleBoxplots <- function(abs, ord, data, color_palette = "Dark2", color_manual = FALSE){

  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  abs <- enquo(abs)
  ord <- enquo(ord)
  
  #Graphic :
  plot <- ggplot(data) +
          #To add points on graph
          geom_quasirandom(aes(y = !!ord, x = !!abs, color = !!abs), size = 1, alpha = 0.5) +
          #To add boxplots
          geom_boxplot(aes(y = !!ord, x = !!abs, color = !!abs, fill = !!abs), outlier.shape = NA, width = 0.4, lwd = 0.6) + 
          #To see points and don't have background color in boxplots
          scale_fill_manual(values = c("transparent", "transparent", "transparent")) +
          #To don't have background color in legend
          guides(fill = "none") +
          #To change background plot
          theme_classic()

  if(color_palette != FALSE){
    #To put color for boxplots with color palette
    plot <- plot + scale_color_brewer(palette = color_palette) 

  }  else {
    #To put color for boxplots with a vector of colors
    plot <- plot + scale_color_manual(values = color_manual)
  }
  
  return(plot)
}

