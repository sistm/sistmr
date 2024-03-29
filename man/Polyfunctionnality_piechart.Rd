% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICS-Polyfunctionnality_piechart.R
\name{Polyfunctionnality_piechart}
\alias{Polyfunctionnality_piechart}
\title{Piecharts for polyfunctionnality of ICS data
 
The complete cytokines must have this type of format : IFN-IL2+MIP1b-TNF+ in CD4}
\usage{
Polyfunctionnality_piechart(
  data_ICS,
  list_cytokines,
  ID_name,
  group_var,
  group_value,
  Timepoint_var,
  Timepoint_value,
  Stim_var,
  Stim_value,
  cellPop_type,
  group_cyto = NULL,
  title_piechart = NULL,
  leg_piechart = TRUE,
  pie_col = NULL,
  ring_col = NULL
)
}
\arguments{
\item{data_ICS}{Dataset to yse for plot.}

\item{list_cytokines}{a vector with the cytokines avalaible in \code{data_ICS} (in order of appareance in boolean cytokines)}

\item{ID_name}{a character to give the colname of subject ID variable in \code{data_ICS}}

\item{group_var}{a character to give the colname of group variable in \code{data_ICS}}

\item{group_value}{a character to select one factor of \code{group_var}}

\item{Timepoint_var}{a character to give the colname of Time point in \code{data_ICS}}

\item{Timepoint_value}{a character to select one factor of \code{Timepoint_var}}

\item{Stim_var}{a character to give the colname of stimulation variable in \code{data_ICS}}

\item{Stim_value}{a character if you want to focus on one stimulation (in case of you have several stimulations in \code{data_ICS})}

\item{cellPop_type}{a character if you have differents populations for your cytokines. It must be present in the column names with cytokines.}

\item{group_cyto}{a vector of numerics if you want to combine two numbers of cytokines (i.e. "c(2,3)" if you want to have a category with "2 and 3" cytokines together). Default is \code{NULL}.}

\item{title_piechart}{a character with your title or \code{NULL} if you don't want a title.}

\item{leg_piechart}{a logical value. If \code{TRUE}, the legend is shown. Default is \code{TRUE}.}

\item{pie_col}{a vector of colors for slices. Default is \code{NULL} then colors of \code{grey.colors} function will be use.}

\item{ring_col}{a vector of colors for rings. Default is \code{NULL} then colors of \code{rainbow} function will be use.}
}
\value{
a \code{pie} object
}
\description{
Piecharts for polyfunctionnality of ICS data
 
The complete cytokines must have this type of format : IFN-IL2+MIP1b-TNF+ in CD4
}
\examples{
Polyfunctionnality_piechart(data_ICS = ICS_data, 
list_cytokines = c("IFN", "IL2", "MIP1b", "TNF"), ID_name = "ID", group_var = "Arm", 
group_value = "Arm1",  Timepoint_var = "TP", Timepoint_value = "TP1", Stim_var = "Stim", 
Stim_value = "BackgroundSubstracted", cellPop_type = "CD4", title_piechart = "Polyfunctionnality", 
leg_piechart = TRUE)
}
\author{
Mélanie Huchon (adapted to Mélany Durand codes)
}
