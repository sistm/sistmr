library(dplyr)
library(reshape2)
library(GSA)
library(illuminaHumanv4.db)
library(ggplot2)

source("normalize.R")


dynamic_genes <- function(data, meta, vec_week, vec_group=NULL,
                  vec_gene,meta_week, meta_group=NULL, meta_sample_ID,
                  meta_participant,meta_sex=NULL, meta_age=NULL, indiv,
                  group_facet=FALSE, convert=FALSE, legend=TRUE) {

  pars <- as.list(match.call()[-1])
  #Mean individus
  if(indiv == 0) {
    #Put data rownames in the first column of the data frame
    df <- tibble::rownames_to_column(as.data.frame(data), "VALUE")
    colnames(df)[1] <- as.character((pars$meta_sample_ID))

    meta[,c(as.character((pars$meta_sample_ID)))] <-
      as.character(meta[,c(as.character((pars$meta_sample_ID)))])
    #Merge expression data and meta data
    merge_data <- merge(meta[,c(as.character((pars$meta_sample_ID)),
                  as.character((pars$meta_participant)),
                  as.character((pars$meta_week)),
                  as.character((pars$meta_group)))],
                  df[,c(1, which(colnames(df) %in% vec_gene))],
                  by=colnames(df)[1])

    #Mean of individus
    mean_indiv_gene <- NULL
    for(visit in sort(vec_week)) {
      #If independant groups
      if(!is.null(pars$meta_group)) {
        for(groupi in sort(vec_group)) {
          if(nrow(merge_data[which(merge_data[,as.character(pars$meta_week)]
            == visit & merge_data[,as.character(pars$meta_group)] == groupi),
            vec_gene]) != 0) mean_indiv_gene <- rbind(mean_indiv_gene,
            c(round(apply(merge_data[which(merge_data[,
            as.character(pars$meta_week)] == visit &
            merge_data[,as.character(pars$meta_group)] == groupi),vec_gene],
            2, mean),
            digits=2), visit, groupi))

        }
      }
      else {
        if(nrow(merge_data[which(merge_data[,as.character(pars$meta_week)] ==
          visit),vec_gene]) != 0) mean_indiv_gene <- rbind(mean_indiv_gene,
          c(round(apply(merge_data[which(merge_data[,
          as.character(pars$meta_week)] == visit),vec_gene],2, mean),
          digits=2), visit))
      }

    }

    mean_indiv_gene <- as.data.frame(mean_indiv_gene)
    colnames(mean_indiv_gene)[length(vec_gene)+1] <-
    as.character(pars$meta_week)
    if(!is.null(pars$meta_group))
    colnames(mean_indiv_gene)[length(vec_gene)+2] <-
      as.character(pars$meta_group)

    mean_indiv_data_plot <- melt(mean_indiv_gene, c(as.character(pars$meta_week
                            ),as.character(pars$meta_group)))

    #Convert ProbeID in illumina ID
    #####TO DO###############
    #To Change Gene column !!
    if(convert==TRUE) {
    probes2ill <- as.list(illuminaHumanv4ARRAYADDRESS)
    ill2symb <- as.list(illuminaHumanv4SYMBOL)
    gene_name <- NULL
    for(i in 1:nrow(mean_indiv_data_plot)) {
      gene_name <-  c(gene_name,as.character(ill2symb[as.character(names(
                    probes2ill[which(probes2ill %in% mean_indiv_data_plot
                    [i,"variable"])]))]))
    }
    mean_indiv_data_plot$Gene <- gene_name
    }
    #################################

    #Normalize
    mean_indiv_data_plot$Norm <- rep(NA,nrow(mean_indiv_data_plot))
    for(gene in unique(mean_indiv_data_plot$variable)) {
      if(!is.null(pars$meta_group)) {
        for(groupi in unique(mean_indiv_data_plot
          [,as.character((pars$meta_group))])) {
          index <- which(mean_indiv_data_plot$variable %in% gene &
                   unique(mean_indiv_data_plot[,as.character((pars$meta_group))]
                   ) %in% groupi)
          mean_indiv_data_plot[index,]$Norm <- normal_distribution(
          mean_indiv_data_plot[index,]$value)
        }
      }
      else {
        index <- which(mean_indiv_data_plot$variable %in% gene)
        mean_indiv_data_plot[index,]$Norm <-
        normal_distribution(mean_indiv_data_plot[index,]$value)
      }
    }

    mean_indiv_data_plot[,as.character(pars$meta_week)] <- as.numeric(
    as.character(mean_indiv_data_plot[,as.character(pars$meta_week)]))

    #Plot
    colnames(mean_indiv_data_plot)[which(colnames
    (mean_indiv_data_plot) %in% "variable")] <- "Gene"
    colnames(mean_indiv_data_plot)[which(colnames
    (mean_indiv_data_plot) %in% as.character(pars$meta_week))] <- "time"
    if(!is.null(pars$meta_group)) colnames(mean_indiv_data_plot)[
    which(colnames(mean_indiv_data_plot) %in% as.character(
    pars$meta_group))] <- "group"

    if(group_facet == TRUE || is.null(pars$meta_group)) {
      p <- ggplot(data=mean_indiv_data_plot,
                  aes(x=time, y = Norm,
                      colour = Gene),na.rm = TRUE)
    }
    if(group_facet == FALSE && !is.null(pars$meta_group)) {
      p <- ggplot(data=mean_indiv_data_plot,
                  aes(x=time, y = Norm, colour = group,
                  group = interaction(group, Gene)),
                  na.rm = TRUE)

    }

    p <- p + geom_point(size=1) +
      geom_line(linetype='dashed') +
      ylab(label = "Gene expression")



    if(!is.null(pars$meta_group)) {
      if(group_facet == TRUE) p <- p + facet_wrap(~ group, ncol = 2)
    }

    p <- p + ggtitle("Dynamic of gene expression") +
         geom_vline(xintercept=unique(mean_indiv_data_plot$time),
         linetype=4, color="#A8A8A8") +
         scale_x_continuous(breaks =
         unique(mean_indiv_data_plot$time)) +
         theme_bw()

    if(legend == FALSE) p <- p + theme(legend.position="none")

    print(p)

    return(mean_indiv_data_plot)

    ##Return plot and data frame
  }

  ##Option with sex and age and group
  else {
    ###TO DO
    #Return plot and data frame
    print("This part is in construction...")
  }


}
