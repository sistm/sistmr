% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ICS_data}
\alias{ICS_data}
\title{Intracellular cytokine staining (ICS) data}
\format{
a data frame with 567 lines and 34 variables:
\describe{
  \item{Random_ID}{Identifiant of participant}
  \item{TP}{Sampling time point}
  \item{Stim}{Type of stimulation ("Non-stimulated", "Stimulated", "BackgroundSubstracted"(Stimulated - non-stimulated))}
  \item{Arm}{Vaccine arm}
  \item{other columns}{represent boolean cytokines}
}
}
\usage{
ICS_data
}
\description{
This dataset contains simulated ICS data at several time point with different vaccine arms.
}
\keyword{datasets}
