#' Piecharts for polyfunctionnality of ICS data
#'  
#' The complete cytokines must have this type of format : IFN-IL2+MIP1b-TNF+ in CD4
#'
#' @param data_ICS Dataset to yse for plot.
#' @param list_cytokines a vector with the cytokines avalaible in \code{data_ICS} (in order of appareance in boolean cytokines)
#' @param ID_name a character to give the colname of subject ID variable in \code{data_ICS}
#' @param group_var a character to give the colname of group variable in \code{data_ICS}
#' @param group_value a character to select one factor of \code{group_var}
#' @param Timepoint_var a character to give the colname of Time point in \code{data_ICS}
#' @param Timepoint_value a character to select one factor of \code{Timepoint_var}
#' @param Stim_var a character to give the colname of stimulation variable in \code{data_ICS}
#' @param Stim_value a character if you want to focus on one stimulation (in case of you have several stimulations in \code{data_ICS})
#' @param cellPop_type a character if you have differents populations for your cytokines. It must be present in the column names with cytokines.
#' @param group_cyto a vector of numerics if you want to combine two numbers of cytokines (i.e. "c(2,3)" if you want to have a category with "2 and 3" cytokines together). Default is \code{NULL}.
#' @param title_piechart a character with your title or \code{NULL} if you don't want a title.
#' @param leg_piechart a logical value. If \code{TRUE}, the legend is shown. Default is \code{TRUE}.
#' @param pie_col a vector of colors for slices. Default is \code{NULL} then colors of \code{grey.colors} function will be use.
#' @param ring_col a vector of colors for rings. Default is \code{NULL} then colors of \code{rainbow} function will be use.
#'
#' @return  a \code{pie} object 
#' 
#' @export
#' 
#' @author Mélanie Huchon (adapted to Mélany Durand codes)
#'
#' @examples
#' Polyfunctionnality_piechart(data_ICS = ICS_data, 
#' list_cytokines = c("IFN", "IL2", "MIP1b", "TNF"), ID_name = "ID", group_var = "Arm", 
#' group_value = "Arm1",  Timepoint_var = "TP", Timepoint_value = "TP1", Stim_var = "Stim", 
#' Stim_value = "BackgroundSubstracted", cellPop_type = "CD4", title_piechart = "Polyfunctionnality", 
#' leg_piechart = TRUE)

Polyfunctionnality_piechart <- function(data_ICS, list_cytokines, ID_name, group_var, group_value, Timepoint_var, Timepoint_value, Stim_var, Stim_value, cellPop_type, group_cyto = NULL, title_piechart = NULL, leg_piechart = TRUE, pie_col = NULL, ring_col = NULL){
  
  #to compute the format for data_ICS()
  data_pie <- dataMngt_piechart(data = data_ICS, 
                            list_cyto = list_cytokines, ID = ID_name, 
                            group = group_var, group_val = group_value, 
                            Timepoint = Timepoint_var, Timepoint_val = Timepoint_value, 
                            Stim = Stim_var, Stim_val = Stim_value, 
                            pop_type = cellPop_type, gp_cyto = group_cyto)
  
  #Focus on group and time point
  
  nb_cyto <- data_pie$nb_cytokine[which(data_pie[,Timepoint_var] == Timepoint_value)]
  index <-  order(data_pie$nb_cyto) # To order according to number of cytokines
  
  nb_cyto <- nb_cyto[index]

  # Extract the values for cytokine activation, handling group_var being NULL
  if (!is.null(group_var)) {
    value_cyto <- data_pie[which(data_pie[, Timepoint_var] == Timepoint_value), group_value][index]
  } else {
    # If no group_var is provided, aggregate the data for all observations
    value_cyto <- data_pie[which(data_pie[, Timepoint_var] == Timepoint_value), "All"]
    value_cyto <- value_cyto[index]
  }
  
  # Create a list to store the cytokine arcs
  list_arcs <- list()
  
  for (cyto in list_cytokines) {
    list_arcs[[cyto]] <- data_pie[which(data_pie[, Timepoint_var] == Timepoint_value), cyto][index]
  }
  
  # Title of the pie chart (use default if not provided)
  if (is.null(title_piechart)) {
    title_piechart <- paste("Piechart for Timepoint", Timepoint_value, "and Stim", Stim_value)
  }
  
  # Create the pie chart using the draw_pie_arcs function
  piechart <- draw_pie_arcs(value = value_cyto, arcs = list_arcs, group = nb_cyto, title = title_piechart, leg = leg_piechart, size = 1, piecolors = pie_col, ringcolors = ring_col)
  
  return(piechart)
}
