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
#' @importFrom stringr str_detect str_count str_flatten
#' @importFrom stats median
#' @export
#' 
#' @keywords internal
#' 
#' @author Mélanie Huchon
#'
#' @examples
#' 
#' d <- dataMngt_piechart(data = ICS_data, list_cyto = c("IFN", "IL2", "MIP1b", "TNF"),
#' ID = "ID", group = "Arm", group_val = "Arm1", Timepoint = "TP", Timepoint_val = "TP1",
#' Stim = "Stim", Stim_val = "BackgroundSubstracted", pop_type = "CD4")

dataMngt_piechart <- function(data, list_cyto, ID, group, group_val = NULL, Timepoint, Timepoint_val = NULL, Stim, Stim_val, pop_type, gp_cyto = NULL){
  
  #########################################################
  # CHECKS
  #########################################################
  
  stopifnot("ID not found : "=any(names(data) == ID))
  stopifnot("group not found : "=any(names(data) == group))
  stopifnot("Timepoint not found : "=any(names(data) == Timepoint))
  stopifnot("Stim not found : "=any(names(data) == Stim))
  
  #_____________________________________________ 1rst function (manage data)
  #########################################################
  # CREATE THE COMBINATIONS OF CYTOKINES 
  #########################################################
  # from list_cyto
  temp <- do.call(cbind.data.frame, lapply(list_cyto, function(x) {
    df <- data.frame(x = c("+", "-"))
    colnames(df) <- x
    return(df)
  }))
  
  combinationCyto <- expand.grid(temp)
  combinationCyto <- data.frame(sapply(colnames(combinationCyto), function(x) combinationCyto[,x] <- paste0(x, combinationCyto[,x])))
  
  CytoList <- do.call(paste, c(combinationCyto, sep="", " in ", pop_type))
  
  #some checks to verify if there are complete cytokines in the data.frame data and warnings if missing complete cytokines 
  if(sum(CytoList %in% colnames(data)) != length(CytoList)){
    warning(paste(CytoList[!CytoList %in% colnames(data)], collapse = ", ", " : Not in data"))
    CytoList <- CytoList[CytoList %in% colnames(data)]
  }
  #_____________________________________________________
  
  #Manage data_ICS
  data_ICS <- data[, c(ID, Timepoint, Stim, group, CytoList)]
  data_ICS <- data_ICS[which(data_ICS[, Timepoint] == Timepoint_val & data_ICS[, group] == group_val), ]
  
  #Create the data.frame 
  table_percent <- as.data.frame(expand.grid(TP = unique(data_ICS[, Timepoint]), combination = CytoList))
  
  #Add if the cytokine is activated in the combination
  table_percent[, as.character(list_cyto)] <- NA
  for(cyto in list_cyto){
    table_percent[, cyto] <- ifelse(str_detect(table_percent$combination, paste0(cyto, "\\+")), 1, 0)
  }
  
  #Count the number of cytokines activated
  table_percent$nb_cytokine <- str_count(table_percent$combination, "\\+")
  
  if(!is.null(gp_cyto)){
    table_percent$nb_cytokine <- ifelse(table_percent$nb_cytokine %in% gp_cyto, str_flatten(gp_cyto, "/"), table_percent$nb_cytokine)
  }
  
  group_vector <- as.character(unique(data_ICS[, group]))
  TP_vector <- unique(data_ICS[, Timepoint])
  
  table_percent[, group_vector] <- NA
  
  #Calculate median for each combination of variables
  for(gp in group_vector){
    for(tp in TP_vector){
      for(cyto in CytoList){
        table_percent[which(table_percent$TP == tp & table_percent$combination == cyto), gp] <- median(data_ICS[which(data_ICS[, group]== gp & data_ICS[, Stim] == Stim_val & data_ICS[, Timepoint] == tp), cyto])
      }
    }
  }
  
  return(table_percent)
}
