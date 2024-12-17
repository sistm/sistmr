#' General FingerPrintplot, for BTM and Chaussabel modules
#'
#' @param res_dgsaseq an object of \code{dgsa_seq} function in \code{dearseq} package.
#' @param data_log2FC  a data.frame with genes in rownames and if (time == NULL) 1 column with Fold Change else several columns are allowed but colnames must contain time's string characters.
#' @param time a string character to filter on time. This time must be present on \code{datalog2FC}'s colnames. Default is \code{NULL}.
#' @param geneset_modules a string character "BTM" or "Chaussabel" to use genesets included or an other gmt object. Default is \code{Chaussabel}.
#'
#' @return a \code{data.frame} object
#' 
#' @keywords internal
#'
table_GS <- function(res_dgsaseq,
                     data_log2FC,
                     time = NULL,
                     geneset_modules) {
  
  try(if(length(time) > 1) 
    stop("time must be a single character string"))
  
  try(if(!is.character(time) & !is.null(time)) 
    stop("time must be a character string, a pattern present in 'data_log2FC' colnames"))
  
  s <- summary(res_dgsaseq, signif_threshold = 1)
  
  res <- data.frame("gene_set" = s$which_signif,
                    "adjusted_pvalues" = s$adj_pval)
  res$gene_set <- as.character(res$gene_set)
  res$description <- NA
  res$median_log2FC <- NA
  res$aggregate <- NA
  
  for(i in 1:nrow(res)){
    
    mod <- res$gene_set[i]
    
    # If mod is a number he is considered as the mod^th geneset of the modules
    if (identical(which(geneset_modules$geneset.names == mod), integer(0))) {
      # get geneset description and cluster
      if (!is.na(as.numeric(mod))) mod <- as.numeric(mod) else stop("check gene_set names")
      res$description[i] <- geneset_modules$geneset.descriptions[mod]
      res$aggregate[i] <- geneset_modules$Cluster[mod]
      rownames(res)[i] <- geneset_modules$geneset.names[mod]
    } else {
      res$description[i] <- geneset_modules$geneset.descriptions[which(geneset_modules$geneset.names == mod)]
      res$aggregate[i] <- geneset_modules$Cluster[which(geneset_modules$geneset.names == mod)]
      rownames(res)[i] <- mod
    }
    # Calculate median log2 Fold Change of the geneset, can be filter by time point
    if (is.null(time)) {
      res$median_log2FC[i] = median(data_log2FC[geneset_modules$genesets[[mod]],],
                                    na.rm = TRUE)
    } else {
      res$median_log2FC[i] = median((data_log2FC[geneset_modules$genesets[[mod]],
                                                 colnames(data_log2FC)[stringr::str_detect(colnames(data_log2FC),
                                                                                           time)]]),
                                    na.rm = TRUE)
    }
  }
  return(res)
}



