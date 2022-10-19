#' Volcano plot function 
#'
#' @param log2fc a variable of the magnitude of change (fold-change) in base log 2 corresponding to the x-axis.
#' @param pValue a variable of statistical significance (p-value) corresponding to the y-axis.
#' @param data a data.frame of differentially expressed results from which the 
#' variable \code{log2fc}, \code{pValue} and \code{geneNames} (if it is used) should be taken. 
#' @param FDR_threshold a threshold of false discovery rate.
#' @param LFC_threshold a threshold of log fold change.
#' @param color a vector of two colors for significant or not significant points.
#' @param geneNames a vector of gene names if you want to put gene tags on the volcano plot. Default is NULL.
#' @param nb_geneTags number of tags for the significant genes if \code{geneNames} 
#' is not NULL. Default is 20 to obtain the tags of the 20 first significant genes.
#' @param logTransformPVal If TRUE, the p-values will have a negative logarithm transformation (base 10). Default is TRUE.
#'
#' @return a \code{ggplot2} object 
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import rlang
#' @importFrom scales trans_new
#' @export
#' 
#' @examples 
#' genes <- paste0("G", 1:500)
#' pval <- runif(500, max = 0.5)
#' log2FC <- runif(500, min = -4, max = 4)
#' 
#' data <- cbind.data.frame(genes, pval, log2FC)
#' 
#' rm(genes, pval, log2FC)
#' volcanoPlot(log2FC, pval, data, geneNames = genes)

volcanoPlot <- function(log2fc, pValue, data, FDR_threshold = 0.05, LFC_threshold = log2(1.5), 
                        color = c("red", "black"), geneNames = NULL, nb_geneTags = 20, logTransformPVal = TRUE){

  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  log2fc_var <- enquo(log2fc)
  pValue_var <- enquo(pValue)
  geneNames_var <- enquo(geneNames)
  
  
  #To obtain the value vector of this following variables which are in 'data'
  log2fc_val <- eval_tidy(log2fc_var, data)
  pValue_val <- eval_tidy(pValue_var, data)
  geneNames_val <- eval_tidy(geneNames_var, data)

  ### checks
  if (!is.numeric(log2fc_val))
    stop("'log2fc' should be numeric")
  if (!is.numeric(pValue_val))
    stop("'pValue' should be numeric")
  # test if pval is between 0 and 1
  if (any(pValue_val < 0 | pValue_val > 1))
    stop("'pValue' should be >= 0 and <= 1") # prevent from being already on log scale
  
  #To add a significant variable in data for differentiate genes which are significant or not on the plot
  data$significant <- "Significant"
  data$significant[pValue_val > FDR_threshold | abs(log2fc_val) < LFC_threshold] <- "Not significant"
  data$significant <- factor(data$significant, 
                             levels = c("Significant", "Not significant"),
                             ordered= TRUE)
  
  # Initial plot
  plot <- ggplot(data = data, aes(x = !!log2fc_var, y = !!pValue_var)) +
          #To add points
          geom_point(data = data, aes_string(color = "significant", fill = "significant"), alpha = 0.5) +
          #To add vertical and horizontal threshold lines
          geom_vline(xintercept = LFC_threshold, lty = 2, colour="black") +
          geom_vline(xintercept = -LFC_threshold, lty = 2, colour="black") +
          geom_hline(yintercept = FDR_threshold, lty = 2, colour="black") +
          #To put the coordinates of the abscissa in log2
          scale_x_continuous(breaks=log2(c(1/50, 1/15, 1/5, 1/1.5, 1.5, 5, 15, 50)), 
                             minor_breaks = NULL,
                             limits = c(-7,7), labels = c("1/50", "1/15", "1/5", "1/1.5", "1.5", "5", "15", "50")) + 
          #To put the log2 scale in the x-axis label
          xlab(paste0(as_name(log2fc_var), " (log2 scale)")) +
          #To put the chosen color and name = "" to remove the title of legend
          scale_color_manual(values = color, name = "") +
          scale_fill_manual(values = color, name = "") +   
          #To change background plot
          theme_bw()
  
  #To obtain the "nb_geneTags" first significant genes and put tags on plot
  if(!is.null(geneNames_val)){
    #To obtain the "nb_geneTags" first significant genes 
    geneTags <- data
    geneTags <- geneTags[which(pValue_val < FDR_threshold),]
    log2fc_gT <- rlang::eval_tidy(log2fc_var, geneTags)
    geneTags <- geneTags[order(abs(log2fc_gT), decreasing = TRUE),]
    geneTags <- geneTags[1:nb_geneTags, ]

    #To add a tags for the significant gene on plot
    data[, "Tags"] <- "No"
    data[which(geneNames_val %in% geneTags$genes_SYMBOL), "Tags"] <- "Yes"
    
    dataToPlot_Tags <- data[which(data$Tags == "Yes"), ]
    
    #plot with tags 
    plot <- plot +
            geom_label_repel(data = dataToPlot_Tags,  
                             aes(x = !!log2fc_var, y = !!pValue_var, label = !!geneNames_var), 
                             min.segment.length = 0, size = 2, alpha = 0.5, seed = 1234, show.legend = FALSE, force = 6, max.overlaps = Inf) +
            geom_label_repel(data =  dataToPlot_Tags, 
                             aes(x = !!log2fc_var, y = !!pValue_var, label = !!geneNames_var), 
                             min.segment.length = 0, size = 2, fill = NA, seed = 1234, show.legend = FALSE, force = 6, max.overlaps = Inf) 
    print(suppressWarnings(plot))
  }
  
  if(logTransformPVal){
    
    #To transform p-values into -log10
    trans_mlog10 <- scales::trans_new(name = "-log10",
                                      transform = function(x){-log10(x + 0.00001)},
                                      inverse = function(u){10^(-u) - 0.00001},
                                      domain = c(0, Inf))
    #To add this transformation on plot
    plot <- plot +
            ylab("FDR p-value (-log10 scale)") +
            scale_y_continuous(trans = trans_mlog10,
                               breaks = c(1,0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005),
                               minor_breaks = c(2:9/10, 2:9/100, 2:9/1000, 2:9/10000), 
                               labels = c(expression("1x10"^{0}), expression("5x10"^{-1}), expression("1x10"^{-1}), 
                                          expression("5x10"^{-2}), expression("1x10"^{-2}), expression("5x10"^{-3}),
                                          expression("1x10"^{-3}), expression("5x10"^{-4}))) 
  } else {
    plot <- plot +
            ylab("FDR p-value") 
  }
  
  return(plot)
  
}