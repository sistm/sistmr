#' Function to manage data for the piecharts of polyfunctionnality of ICS data  
#' 
#' The boolean cytokines must have this type of format : IFN-IL2-MIP1b-TNF+ in CD4
#' 
#' @param data 
#' @param list_cyto 
#' @param ID a character to define the colnames of ID in data_ICS
#' @param group a character to define the colnames of group in data_ICS (Arm)
#' @param group_val default is NULL : no filter
#' @param Timepoint a character to define the colnames of Time point in data_ICS 
#' @param Timepoint_val default is NULL : no filter
#' @param Stim a character to define the colnames of stimulation in data_ICS
#' @param Stim_val 
#' @param pop_type a character to define the type of population (CD4 or CD8) (if your cytokines are for example : "IFN+IL2+MIP1B+ in CD4" in data)
#' @param gp_cyto vector with the number of cytokines that you want to group (ex : c(2,3))
#'
#' @return a \code{data.frame} object
#' 
#' @importFrom stringr str_detect str_count str_flatten str_which str_replace_all
#' @importFrom stats median
#' 
#' @keywords internal
#' 
#' @author Mélanie Huchon
#'
#' @examples
#' 
#' \dontrun{d <- dataMngt_piechart(data = ICS_data, list_cyto = c("IFN", "IL2", "MIP1b", "TNF"),
#' ID = "ID", group = "Arm", group_val = "Arm1", Timepoint = "TP", Timepoint_val = "TP1",
#' Stim = "Stim", Stim_val = "BackgroundSubstracted", pop_type = "CD4")}

dataMngt_piechart <- function(data, list_cyto, ID, group, group_val = NULL, Timepoint, Timepoint_val = NULL, Stim, Stim_val, pop_type, gp_cyto = NULL){
  
  #########################################################
  # CHECKS
  #########################################################
  
  # Check that "ID" and "Stim" columns exist in the dataset
  stopifnot("ID not found : " = any(names(data) == ID))
  stopifnot("Stim not found : " = any(names(data) == Stim))
  
  # If the data is a tibble, convert it to a data.frame for compatibility
  if(is.tbl(data)){
    data <- as.data.frame(data)
  }
  
  #_____________________________________________ 1rst function (manage data)
  #########################################################
  # CREATE THE COMBINATIONS OF CYTOKINES 
  #########################################################
  # Create all possible cytokine combinations (active or not) from list_cyto
  temp <- do.call(cbind.data.frame, lapply(list_cyto, function(x) {
    df <- data.frame(x = c("+", "-"))  # For each cytokine, create two states: "+" (active) or "-" (inactive)
    colnames(df) <- x
    return(df)
  }))
  
  # Create a table of all possible cytokine combinations
  combinationCyto <- expand.grid(temp)
  # Format combinations to include cytokine names with their status (+ or -)
  combinationCyto <- data.frame(sapply(colnames(combinationCyto), function(x) combinationCyto[,x] <- paste0(x, combinationCyto[,x])))
  
  # Combine the cytokine combinations into a single string for each row
  combinationCyto_list <- do.call(paste, c(combinationCyto, sep=""))
  
  # Use the find_column_names function to retrieve column names matching these combinations
  CytoList <- find_column_names(column_names = colnames(data), population = pop_type, cytokine_combinations = combinationCyto_list)
  
  # Check if all cytokine combinations are present in the dataset columns
  if(sum(CytoList %in% colnames(data)) != length(CytoList)){
    warning(paste(CytoList[!CytoList %in% colnames(data)], collapse = ", ", " : Not in data"))
  }
  #_____________________________________________________
  
  # Manage data_ICS
  # Select relevant columns (ID, Timepoint, Stim, and cytokine combinations)
  selected_columns <- c(ID, Timepoint, Stim, CytoList)
  if (!is.null(group)) {
    selected_columns <- c(selected_columns, group)  # Add group column if not NULL
  }
  
  # Create a subset of data containing only the selected columns
  data_ICS <- data[, selected_columns]
  # Filter rows where the Timepoint value matches Timepoint_val
  data_ICS <- data_ICS[which(data_ICS[, Timepoint] == Timepoint_val), ]
  
  # If a group is specified, filter rows where the group value matches group_val
  if(!is.null(group)){
    data_ICS <- data_ICS[which(data_ICS[, group] == group_val), ]
  }
  
  # Create a new data.frame with all combinations of Timepoint and cytokines
  table_percent <- as.data.frame(expand.grid(TP = unique(data_ICS[, Timepoint]), combination = CytoList))
  colnames(table_percent)[1] <- Timepoint

  # Add columns to indicate if each cytokine is active (+) or not (-) in the combination
  table_percent[, as.character(list_cyto)] <- NA
  for(cyto in list_cyto){
    table_percent[, cyto] <- ifelse(str_detect(table_percent$combination, paste0(cyto, "\\+")), 1, 0)
  }
  
  # Count the total number of active cytokines for each combination
  table_percent$nb_cytokine <- rowSums(table_percent[, as.character(list_cyto)], na.rm = TRUE)
  
  # If gp_cyto is specified, adjust the nb_cytokine column according to the groups defined in gp_cyto
  if(!is.null(gp_cyto)){
    table_percent$nb_cytokine <- ifelse(table_percent$nb_cytokine %in% gp_cyto, str_flatten(gp_cyto, "/"), table_percent$nb_cytokine)
  }
  
  # If group is not NULL, create a group vector
  if (!is.null(group)) {
    group_vector <- as.character(unique(data_ICS[, group]))
  } else {
    group_vector <- "All"  # Otherwise, use a default group called "All"
  }
  
  # Create a vector of unique Timepoints
  TP_vector <- unique(data_ICS[, Timepoint])
  
  # Add a column for each group, initialized to NA
  table_percent[, group_vector] <- NA
  
  # Calculate the median for each combination of variables (Timepoint, group, cytokine)
  for(gp in group_vector){
    for(tp in TP_vector){
      for(cyto in CytoList){
        if (gp == "All") {
          # If group is NULL, calculate the median for all data matching Timepoint and Stim
          table_percent[which(table_percent[, Timepoint] == tp & table_percent$combination == cyto), gp] <- median(data_ICS[which(data_ICS[, Stim] == Stim_val & data_ICS[, Timepoint] == tp), cyto], na.rm = TRUE)
        } else {
          # Otherwise, calculate the median by group
          table_percent[which(table_percent[, Timepoint] == tp & table_percent$combination == cyto), gp] <- median(data_ICS[which(data_ICS[, group] == gp & data_ICS[, Stim] == Stim_val & data_ICS[, Timepoint] == tp), cyto], na.rm = TRUE)
        }
      }
    }
  }
  
  return(table_percent)  # Return the final table with calculated percentages/medians
}


find_column_names <- function(column_names, population, cytokine_combinations) {
  results <- vector()  # Initialize a vector to store results
  
  # Check for columns containing the population name
  population_present <- stringr::str_which(column_names, population)
  
  cytokine_match_found <- FALSE  # Variable to indicate if combinations are found
  for (cytokine in cytokine_combinations) {
    # Format cytokine combination for use with stringr functions
    cytokine_combination_str <- str_replace_all(cytokine, "\\+", "\\\\+")
    
    # Find column names that match the cytokine combination
    index_colnames_match <- str_which(column_names, cytokine_combination_str)
    
    # Append the found column names to the result vector
    cytokine_match_found <- c(index_colnames_match, cytokine_match_found)
  }
  
  # Get the intersection of columns containing the population and matching cytokines
  results <- column_names[intersect(population_present, cytokine_match_found)]
  
  return(results)  # Return the found column names
}