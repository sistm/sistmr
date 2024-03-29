% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colorDendro.R
\name{colorDendro}
\alias{colorDendro}
\title{Color the dendrogram labels}
\usage{
colorDendro(hclust, labels = hclust$labels, lab.col = NULL, ...)
}
\arguments{
\item{hclust}{a \code{hclust} object.}

\item{labels}{a character vector of labels for the leaves of the tree. By default, the labels in the \code{hclust} object are used if specified.}

\item{lab.col}{a color vector for the labels.}

\item{...}{further arguments passed to plot function}
}
\value{
a dendrogram
}
\description{
Assign colors to the labels of a \code{hclust} dendrogram according to a factor
}
\examples{
#Load data example
data(iris)
#Run hclust on this data
hc <- hclust(dist(iris[, -5]))
#Color vector 
color <- c("#440154FF", "#21908CFF", "#FDE725FF")
#Use function to color labels
colorDendro(hc, labels = iris[,5], lab.col = color, main = "Dendrogram of Iris data")
}
\author{
Mélanie Huchon
}
