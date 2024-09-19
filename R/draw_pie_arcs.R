#' Draw pie charts with arcs
#' 
#' This function is to do plots for polyfunctionnality of ICS data
#'
#' @param value 
#' @param arcs 
#' @param group 
#' @param title 
#' @param leg 
#' @param size 
#' @param pie_col 
#' @param ring_col 
#'
#' @return a \code{pie} object
#'  
#' @importFrom grDevices grey.colors rainbow recordPlot
#' @importFrom graphics pie legend plot.window polygon 
#' 
#' @keywords internal
#' 
#' @author Claire Bauduin (adapted to Binbin Xu codes) - Last modifications by Mélanie Huchon
#'

draw_pie_arcs <- function(value, arcs, group, title, leg = F, size, piecolors = NULL, ringcolors = NULL) {
  
  if(is.null(piecolors)){
    piecolors <- grey.colors(n = length(unique(group)), start = 0, end = 0.95, rev = T, gamma = 1)# colors for sclices
  }
  if(is.null(ringcolors)){
    ringcolors <- rainbow(n = length(arcs), end  = 0.5)# color for rings
  }
 

  ringval <- value / sum(value)
  anglval <- (pi / 2) + c(0, cumsum(ringval) * (2 * pi))
  
  xlim <- ylim <- c(-3, 3)
  plot.window(xlim, ylim, "", asp = 1)
  par(mar = c(0,0,0,0))
  pie(value, radius = 0.5, labels = {''}, col = piecolors[as.numeric(factor(group))], border = NA, init.angle = 90, main = title, cex.main = size)
  # ajout de la l?gende avec les modalit?s : 1/2/3 cytokines et IFN/IL2/TNF
  if (leg)
    legend_pie <- legend("topright", c(paste(unique(group), "cytokines"), names(arcs)), cex = size, fill = c(piecolors,ringcolors))
  
  for (i in 1:length(arcs)) {
    radius0 = 1.7 + (i - 1) * 0.13
    radius1 = 1.7 + (i - 1) * 0.13 + 0.1
    
    cl = ringcolors[i]
    
    for (j in 1:length(value)) {
      if (arcs[[i]][j] == 1) {
        theta0 = anglval[j]
        theta1 = anglval[j + 1]
        
        a1 <- seq(theta0,  theta0,  length = 100)
        r1 <- seq(radius0, radius1, length = 100)
        a2 <- seq(theta0,  theta1,  length = 100)
        r2 <- seq(radius1, radius1, length = 100)
        a3 <- seq(theta1,  theta1,  length = 100)
        r3 <- seq(radius1, radius0, length = 100)
        a4 <- seq(theta1,  theta0,  length = 100)
        r4 <- seq(radius0, radius0, length = 100)
        
        xy <- pol2cart(c(a1, a2, a3, a4), c(r1, r2, r3, r4))
        
        #par(plt = gridPLT(), new = TRUE, mar = c(0,0,0,0))
        plot.window(xlim, ylim, "", asp = 1)
        polygon(xy[, 1], y = xy[, 2], col = cl, border = NA)
      }
    }
  }
  plot <- recordPlot()
  
  return(plot)
}
