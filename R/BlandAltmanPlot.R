#' Bland-Altman plot function
#' 
#' @param var1 a vector of numerics for the 1rst group to be compared.
#' @param var2 a vector of numerics for the 2nd group to be compared.
#' @param with_gradient a logical indicating if you have a lot of measures, use \code{with_gradient=TRUE}
#' to have gradient scale and not points. Default value is FALSE.
#' @param line_color a vector of color for the three lines : average difference and
#'  upper and lower limits of the confidence interval for the average difference.
#' @param extremum_pctg a logical indicating if you want to add the percentage of 
#' points outside the confidence interval for the upper and lower limits. Default is TRUE.
#' 
#' @return a \code{ggplot2} object
#' @import BlandAltmanLeh
#' @import ggplot2
#' @export


BlandAltmanPlot <- function(var1, var2, with_gradient = FALSE, line_color = c("blue", "lightblue"), extremum_pctg = TRUE){

  #To calculate statistics for Bland-Altman plot with default parameters
  ba.stats <- bland.altman.stats(group1 = var1, group2 = var2)
  
  #Create a dataframe for plotting with ggplot2
  data_plot <- cbind.data.frame(means = ba.stats$means, diffs = ba.stats$diffs)
  
  #Initial plot
  plot <- ggplot(data_plot, aes_string(x = "means", y = "diffs")) + 
          theme_classic()
  
  #To add percentages of dots outside the limits (upper/lower)
  if(extremum_pctg){
    #To compute these percentages
    nb_inf <- 0
    nb_sup <- 0
    for (i in 1:length(ba.stats$diffs)){
      if(ba.stats$diffs[i] > ba.stats[["lines"]][["upper.limit"]]){
        nb_sup <- nb_sup + 1
      }else{
        nb_sup <- nb_sup
      }
      if(ba.stats$diffs[i] < ba.stats[["lines"]][["lower.limit"]]){
        nb_inf <- nb_inf + 1
      }else{
        nb_inf <- nb_inf
      }
    }
    pctg_sup <- round((nb_sup/length(ba.stats$diffs))*100, 2)
    pctg_inf <- round((nb_inf/length(ba.stats$diffs))*100, 2)
    
    #To put on the plot
    plot <- plot +
            annotate("text", x = max(data_plot$means)/2, y = max(data_plot$diffs), label = paste0(pctg_sup," % sup")) +
            annotate("text", x = max(data_plot$means)/2, y = min(data_plot$diffs), label = paste0(pctg_inf," % inf")) 
  }
  
  #To put a gradient scale when there are a lot of points
  if(with_gradient){
    plot <- plot + 
            geom_hex(bins = 100) +
            scale_fill_gradient(low = "white", high = "black", trans = "log10")
  } else {
    #To add points 
    plot <- plot +
            geom_point() 
  }
  
  #To add lines and labels 
  plot <- plot +
          #To add lines for the extremum limits and median value
          geom_hline(yintercept = ba.stats$lines[1], linetype = 2, color = line_color[2], size = 1) +
          geom_hline(yintercept = ba.stats$lines[2], linetype = 5, color = line_color[1], size = 1) +
          geom_hline(yintercept = ba.stats$lines[3], linetype = 2, color = line_color[2], size = 1) +
          #To put axis labels 
          xlab("Mean") + 
          ylab("Difference")
  
  return(plot)
  
}