% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcanoPlot.R
\name{volcanoPlot}
\alias{volcanoPlot}
\title{Volcano plot function}
\usage{
volcanoPlot(
  log2fc,
  pValue,
  data,
  FDR_threshold = 0.05,
  LFC_threshold = log2(1.5),
  color = c("red", "black"),
  geneNames = NULL,
  nb_geneTags = 20,
  logTransformPVal = TRUE
)
}
\arguments{
\item{log2fc}{a variable of the magnitude of change (fold-change) in base log 2 corresponding to the x-axis.}

\item{pValue}{a variable of statistical significance (p-value) corresponding to the y-axis.}

\item{data}{a data.frame of differentially expressed results from which the 
variable \code{log2fc}, \code{pValue} and \code{geneNames} (if it is used) should be taken.}

\item{FDR_threshold}{a threshold of false discovery rate.}

\item{LFC_threshold}{a threshold of log fold change.}

\item{color}{a vector of two colors for significant or not significant points.}

\item{geneNames}{a vector of gene names if you want to put gene tags on the volcano plot. Default is NULL.}

\item{nb_geneTags}{number of tags for the significant genes if \code{geneNames} 
is not NULL. Default is 20 to obtain the tags of the 20 first significant genes.}

\item{logTransformPVal}{If TRUE, the p-values will have a negative logarithm transformation (base 10). Default is TRUE.}
}
\value{
a \code{ggplot2} object
}
\description{
Volcano plot function
}
\examples{
genes <- paste0("G", 1:500)
pval <- runif(500, max = 0.5)
log2FC <- runif(500, min = -4, max = 4)

data <- cbind.data.frame(genes, pval, log2FC)

rm(genes, pval, log2FC)
volcanoPlot(log2FC, pval, data, geneNames = genes)
}
