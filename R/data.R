#' Intracellular cytokine staining (ICS) data
#'
#' This dataset contains simulated ICS data at several time point with different vaccine arms.
#'
#' @format a data frame with 567 lines and 34 variables:
#' \describe{
#'   \item{Random_ID}{Identifiant of participant}
#'   \item{TP}{Sampling time point}
#'   \item{Stim}{Type of stimulation ("Non-stimulated", "Stimulated", "BackgroundSubstracted"(Stimulated - non-stimulated))}
#'   \item{Arm}{Vaccine arm}
#'   \item{other columns}{represent boolean cytokines}
#' }
"ICS_data"