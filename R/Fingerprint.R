#' Function for Finger print plot on Gene set analysis 
#' 
#' The aim of this function is to create a visual representation of gene set enrichment results from differential gene expression analysis.
#' 
#' @param res_dgsaseq an object of \code{dgsa_seq} function in \code{dearseq} package.
#' @param datalog2FC a data.frame with genes in rownames and if (time == NULL) 1 column with Fold Change else several columns are allowed but colnames must contain time's string characters.
#' @param genesets a string character "BTM" or "Chaussabel" to use genesets included or an other gmt object. Default is \code{Chaussabel}.
#' @param category a string character or vector. Only needed if \code{(genesets == "BTM")} : "immune", "signaling", "TF targets", "biological process", "molecular function", "location", "all".
#' @param time a string character to filter on time. This time must be present on \code{datalog2FC}'s colnames. Default is \code{NULL}.
#' @param alpha a numeric corresponding to type 1 error threshold. Default value is \code{0.05}.
#' @param title a string character for plot's title. Default is \code{NULL}.
#' @param max_size a numeric corresponding to the maximum size of point. Default is \code{8}.  
#' @param size_y_lab labsize for non-interactive plot, can be useful to set to 0 to combine several plots
#' @param order_by_pvals a logical indicating to order genesets per pvalue left to right. Default is \code{FALSE}.
#' @param interactive.plot a logical indicating whether the plot will be interactive, recommended if order_by_pvals is \code{TRUE}. Default is \code{FALSE}.
#'
#' @return a \code{ggplot} object
#' 
#' @importFrom stats setNames
#' @importFrom ggiraph geom_point_interactive opts_tooltip girafe
#' @importFrom scales rescale
#' @importFrom reshape2 melt
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stringr str_replace
#' 
#' @author Mélanie Huchon & Quentin Laval
#' 
#' @export
#'
#' @examples
#' 

Fingerprint <- function(res_dgsaseq,
                        datalog2FC,
                        genesets = 'Chaussabel',
                        category = 'immune',
                        time = NULL,
                        alpha = 0.05,
                        title = NULL,
                        max_size = 8,
                        size_y_lab = 9, 
                        order_by_pvals = FALSE,
                        interactive.plot = FALSE){

  stopifnot("Please select Chaussabel / BTM or put a personnal geneset"=
              is.character(genesets) & genesets %in% c("Chaussabel", "BTM") |
              is.list(genesets))
  
  if (all(!genesets %in% c("Chaussabel", "BTM"))) {
    stopifnot(
      "Your personal geneset must contain the following list names: genesets, geneset.names, geneset.descriptions and Cluster" =
        all(c("genesets", "geneset.names", "geneset.descriptions", "Cluster") %in% names(genesets))
    )
  }
  stopifnot("'order_by_pvals' must be TRUE or FALSE"=order_by_pvals %in% c(TRUE, FALSE))
  stopifnot("'interactive.plot' must be TRUE or FALSE"=interactive.plot %in% c(TRUE, FALSE))
  
  if (all(is.character(genesets) & genesets == "Chaussabel")) {
    
    load(file = "data/gmtChaussabel_sistmr.RData")
    gmt <- gmtChaussabel
    res_geneset <- table_GS(res_dgsaseq = res_dgsaseq,
                            time = time,
                            data_log2FC = datalog2FC,
                            geneset_modules = gmt)
    # Remove TBD modules
    id_TBD <- which(gmt$geneset.descriptions == 'TBD')
    gmt$geneset.descriptions <- gmt$geneset.descriptions[- id_TBD]
    gmt$geneset.names <- gmt$geneset.names[- id_TBD]
    if (!identical(which(res_geneset$description == 'TBD'), integer(0))){
      
      res_geneset <- res_geneset[-which(res_geneset$description == 'TBD'),]
      
    }
  } else if (all(is.character(genesets) & genesets == "BTM")) {
    
    load(file = "data/gmtBTM_sistmr.RData")
    gmt <- gmtBTM
    res_geneset <- table_GS(res_dgsaseq = res_dgsaseq,
                            time = time,
                            data_log2FC = datalog2FC,
                            geneset_modules = gmt)
    
    BTM_categories <- c("immune", "signaling", "TF targets", "biological process",
                        "molecular function", "location", "all")
    
    if (!all(category %in% BTM_categories)) {
      stop(paste("Please select one or more categories : ",
                 paste(BTM_categories, collapse = ", ")))
    }
    
    # Filter gmt according to category 
    if (length(which(category == "all")) > 0) {
      sub_mods <- BTM_categories[-7]
    } else {
      
      sub_mods <- category
      match_category <- which(gmt$Cluster %in% category)
      gmt <- lapply(setNames(names(gmt), names(gmt)), function(b)
        gmt[[b]] <- gmt[[b]][match_category])
    }
    # Subset modules via category
    res_geneset <- res_geneset[which(res_geneset$aggregate %in% sub_mods),]
  } else {
    
    gmt <- genesets
    res_geneset <- table_GS(res_dgsaseq = res_dgsaseq,
                            time = time,
                            data_log2FC = datalog2FC,
                            geneset_modules = gmt)
    if (!is.null(category)){
      match_category <- which(gmt$Cluster %in% category)
      gmt <- lapply(setNames(names(gmt), names(gmt)), function(b)
        gmt[[b]] <- gmt[[b]][match_category])
      # Subset modules via category
      res_geneset <- res_geneset[which(res_geneset$aggregate %in% category),]
    } else {
      # Remove TBD modules
      id_TBD <- which(gmt$geneset.descriptions == 'TBD')
      gmt$geneset.descriptions <- gmt$geneset.descriptions[- id_TBD]
      gmt$geneset.names <- gmt$geneset.names[- id_TBD]
      if (!identical(which(res_geneset$description == 'TBD'), integer(0))){
        
        res_geneset <- res_geneset[-which(res_geneset$description == 'TBD'),]
        
      }
    }
  }
  
  # Name of functionnal clusters and right rows names
  cluster_names <- unique(gmt$geneset.descriptions)
  
  # Prepare data
  missing_modules <- which(is.na(match(gmt$geneset.names,
                                       rownames(res_geneset))))
  if (!identical(missing_modules, integer(0))) {
    
    res_geneset <- rbind(res_geneset,
                         data.frame(gene_set = missing_modules,
                                    adjusted_pvalues = NA,
                                    description = gmt$geneset.descriptions[missing_modules],
                                    median_log2FC = NA,
                                    aggregate = res_geneset$aggregate[1],
                                    row.names = gmt$geneset.names[missing_modules]))
  }
  res_geneset$gene_set <- rownames(res_geneset)
  Group_plot <- res_geneset
  
  # Create a variable "Stat" which combine the pvalue and FC informations
  Group_plot$Stat <- Group_plot$adjusted_pvalues
  Group_plot <- Group_plot[, "Stat", drop = FALSE]
  Group_plot <- Group_plot[as.character(rownames(res_geneset)),
                           1, drop = FALSE]
  
  Group_plot <- as.data.frame(Group_plot)
  Group_plot[is.na(Group_plot) == TRUE] <- 1
  
  # Create new grid with all filtered cluster
  grid_mat <- matrix(nrow = length(cluster_names),
                     ncol = max(table(res_geneset$description)))
  rownames(grid_mat) = sort(cluster_names)
  colnames(grid_mat) = seq_len(max(table(res_geneset$description)))
  
  # Extract all x and y coordinates and counter creation for y coord 
  num_desc <- as.numeric(factor(gmt$geneset.descriptions,
                                levels = sort(cluster_names),
                                labels = sort(cluster_names)))
  counter <- integer(length(cluster_names))
  
  gmt$Cluster_position <- 0
  for (p in seq_along(gmt$geneset.descriptions)) {
    
    val <- num_desc[p]
    counter[val] <- counter[val] + 1
    gmt$Cluster_position[p] <- paste(val, counter[val], sep = ".")
    
  }
  match_id <- match(rownames(res_geneset), gmt$geneset.names) 
  res_geneset$Cluster_position <- gmt$Cluster_position[match_id]
  
  
  # Order Cluster_position by adj_pvalues rank in each cluster
  if (order_by_pvals) {
    
    res_geneset <- res_geneset %>%
      rownames_to_column(var = "rownames") %>%  # Convert rownames to a column
      group_by(description) %>%
      mutate(rank = rank(adjusted_pvalues, ties.method = "first")) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(
        Cluster_position = str_replace(
          Cluster_position,
          "(?<=\\.)[0-9]+", 
          as.character(rank)
        )
      ) %>%
      ungroup() %>%
      select(-rank) %>%
      column_to_rownames(var = "rownames")  # Restore rownames
    
  }
  
  x_coords <- res_geneset[rownames(Group_plot), "description"]
  y_coords <- sapply(strsplit(res_geneset[rownames(Group_plot),
                                          "Cluster_position"],
                              "\\."), function(x) x[[2]])
  
  # Extract assign values
  values <- Group_plot$Stat
  
  # Convert x and y coordinates into matrix indices
  x_indices <- match(x_coords, rownames(grid_mat))
  y_indices <- as.numeric(y_coords)
  
  # Update grid_mat with the corresponding values
  grid_mat[cbind(x_indices, y_indices)] <- values
  
  melted_df <- melt(grid_mat, id.var = "row.names")
  colnames(melted_df) <- c("Aggregate", "Sub_aggregate", "AdjPval")
  melted_df$AdjPval <- abs(melted_df$AdjPval)
  
  res_geneset$Cluster_join <- paste0(res_geneset$description,
                                     ".",
                                     lapply(strsplit(res_geneset$Cluster_position,"\\."),
                                            function(x) x[[2]]) )
  
  res_geneset$gene_setID <- res_geneset$gene_set
  
  melted_df <- left_join(melted_df %>%
                           mutate(Cluster_join = paste0(Aggregate, ".", Sub_aggregate)),
                         res_geneset[, c("Cluster_join", "gene_setID")],
                         by = "Cluster_join") %>%
    select(-Cluster_join)
  
  melted_df <- left_join(melted_df,
                         res_geneset[, c("gene_set", "adjusted_pvalues",
                                         "median_log2FC")],
                         join_by(gene_setID == gene_set,
                                 AdjPval == adjusted_pvalues))

  # Tile fill color
  color <- rep("white", nrow(melted_df))
  color[is.na(melted_df$AdjPval)] <- "#E6E6E6"
  recurrence <- length(cluster_names)
  
  # To create an empty white tile when an existing module isn't DE 
  # "#E6E6E6" equals to na.value for the plot
  gray_id <- which(color == "#E6E6E6")
  
  # Shifted indices by recurrence
  white_id <- gray_id + recurrence
  
  # Make sure offset indices are within the limits of the vector
  valid_id <- white_id <= length(color)
  
  # Apply changes only on valid_id
  color[gray_id[valid_id]] <- ifelse(color[white_id[valid_id]] == "white",
                                     "white",
                                     color[gray_id[valid_id]])
  
  melted_df$AdjPval[gray_id[valid_id]] <- ifelse(color[white_id[valid_id]] == "white",
                                                 1,
                                                 melted_df$AdjPval[gray_id[valid_id]])
  
  melted_df$signif <- melted_df$AdjPval <= alpha

  if (interactive.plot) {
    # Remove tooltip for NA
    id <- which(is.na(melted_df$gene_setID))
    melted_df$gene_setID[id] <- ""
    
    melted_df$tooltip_text <- paste(
      melted_df$gene_setID,
      sprintf("FC: %.2f", 2^melted_df$median_log2FC),
      sprintf("FDR: %.2e", melted_df$AdjPval),
      sep = "\n"
    )
    
    plot <- ggplot(melted_df, aes(Aggregate, as.factor(Sub_aggregate))) + 
      geom_tile(color = "grey85", linewidth = 0.2, fill = color) +
      geom_point_interactive(aes(size = AdjPval, colour = median_log2FC, shape = signif,
                     tooltip = gene_setID), show.legend = TRUE) +
      scale_size(name = "FDR",
                 breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 1),
                 range = c(max_size, 0),
                 transform = "sqrt",
                 labels = c(0.001, 0.01, 0.05, 0.1, 0.2, 1),
                 limits = c(min(min(melted_df$AdjPval, na.rm = TRUE), 10^-3), 1)) +
      scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19), na.translate = FALSE, name = "Significant FDR") +
      guides(size = guide_legend(reverse=F, override.aes = list(shape = c(rep(15, 3), 
                                                                          rep(19,3)))),
             shape = "none") +
      scale_color_gradientn(colours = c("#002C63", "darkblue", "blue", "white", "red", "darkred", "#660000"),
                            values = rescale(x = log2(c(1/3, 2/3, 1, 1.5, 3))),
                            limits = c(-log2(3.1), log2(3.1)),
                            breaks = log2(c(1/3, 1/2, 2/3, 1, 1.5, 2, 3)),
                            labels = c("1/3", "1/2", "2/3", "1", "1.5", "2", "3"),
                            name = "Fold-Change") +
      ylab("") +
      xlab("") +
      labs(x = "", y = "", title = title) + 
      theme_light() + 
      theme(panel.grid.minor = element_line(colour = "black", 
                                            linewidth = 0.32)) + 
      coord_flip() + 
      scale_x_discrete(limits = rev(levels(melted_df$Aggregate))) + 
      theme(panel.border = element_rect(color = "black", linewidth = 0.18), 
            axis.text.x = element_blank(), 
            axis.text.y = element_text(colour = "black", size = 6,
                                       angle = 0, hjust = 1,
                                       vjust = 0.5, face = "plain"),
            title = element_text(colour = "black" , size = 6),
            legend.text = element_text(colour = "black", size = 6))
    plot <- girafe(ggobj = plot, options = list(
      opts_tooltip(offx = 40,
                   offy = 40,
                   css = "font-size: 45px;")),
      pointsize = 12)
    
  } else {
    
    plot <- ggplot(melted_df, aes(Aggregate, as.factor(Sub_aggregate))) + 
      geom_tile(color = "grey85", linewidth = 0.2, fill = color) +
      geom_point(aes(size = AdjPval, colour = median_log2FC, shape = signif), show.legend = TRUE) +
      scale_size(name = "FDR",
                 breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 1),
                 range = c(max_size, 0),
                 transform = "sqrt",
                 labels = c(0.001, 0.01, 0.05, 0.1, 0.2, 1),
                 limits = c(min(min(melted_df$AdjPval, na.rm = TRUE), 10^-3), 1)) +
      scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19), na.translate = FALSE, name = "Significant FDR") +
      guides(size = guide_legend(reverse=F, override.aes = list(shape = c(rep(15, 3), 
                                                                          rep(19,3)))),
             shape = "none") +
      scale_color_gradientn(colours = c("#002C63", "darkblue", "blue", "white", "red", "darkred", "#660000"),
                            values = rescale(x = log2(c(1/3, 2/3, 1, 1.5, 3))),
                            limits = c(-log2(3.1), log2(3.1)),
                            breaks = log2(c(1/3, 1/2, 2/3, 1, 1.5, 2, 3)),
                            labels = c("1/3", "1/2", "2/3", "1", "1.5", "2", "3"),
                            name = "Fold-Change") +
      ylab("") +
      xlab("") +
      labs(x = "", y = "", title = title) + 
      theme_light() + 
      theme(panel.grid.minor = element_line(colour = "black", 
                                            linewidth = 0.32)) + 
      coord_flip() + 
      scale_x_discrete(limits = rev(levels(melted_df$Aggregate))) + 
      theme(panel.border = element_rect(color = "black", linewidth = 0.18), 
            axis.text.x = element_blank(), 
            axis.text.y = element_text(colour = "black" , size = size_y_lab,
                                       angle = 0, hjust = 1,
                                       vjust = 0.5, face = "plain"))
    
  }
  
  return(plot)
}
