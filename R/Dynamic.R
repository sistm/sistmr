rm(list=ls())
library(dplyr)
library(reshape2)
library(GSA)
library(illuminaHumanv4.db)
library(ggplot2)

source("normalize.R")


dynamic_genes <- function(data, meta, vec_week, vec_group=NULL, vec_gene,meta_week, meta_group=NULL, meta_sample_ID, meta_participant,meta_sex=NULL, meta_age=NULL, indiv) {
  pars <- as.list(match.call()[-1])
  a <- 2
  #Mean individus
  if(indiv == 0) {
    #Put data rownames in the first column of the data frame
    df <- tibble::rownames_to_column(as.data.frame(data), "VALUE")
    colnames(df)[1] <- as.character((pars$meta_sample_ID))

    meta[,c(as.character((pars$meta_sample_ID)))] <- as.character(meta[,c(as.character((pars$meta_sample_ID)))])
    #Merge expression data and meta data
    merge.data <- merge(meta[,c(as.character((pars$meta_sample_ID)),as.character((pars$meta_participant)),as.character((pars$meta_week)),as.character((pars$meta_group)))],
          df[,c(1, which(colnames(df) %in% vec_gene))], by=colnames(df)[1])

    #Mean of individus
    mean.indiv.gene <- NULL
    for(visit in sort(vec_week)) {
      #If independant groups
      if(!is.null(meta_group)) {
        for(group in sort(vec_group)) {
          if(nrow(merge.data[which(merge.data[,as.character(pars$meta_week)] == visit &
            merge.data[,as.character(pars$meta_group)] == group),vec_gene]) != 0) mean.indiv.gene <- rbind(mean.indiv.gene,
              c(round(apply(merge.data[which(merge.data[,as.character(pars$meta_week)] == visit &
                merge.data[,as.character(pars$meta_group)] == group),vec_gene],2, mean), digits=2), visit, group))

        }
      }
      else {
        if(nrow(merge.data[which(merge.data[,as.character(pars$meta_week)] == visit),vec_gene]) != 0) mean.indiv.gene <- rbind(mean.indiv.gene,
           c(round(apply(merge.data[which(merge.data[,as.character(pars$meta_week)] == visit),vec_gene],2, mean), digits=2), visit))
      }

    }

    mean.indiv.gene <- as.data.frame(mean.indiv.gene)
    colnames(mean.indiv.gene)[4] <- as.character(pars$meta_week)
    if(!is.null(meta_group)) colnames(mean.indiv.gene)[5] <- as.character(pars$meta_group)

    mean.indiv.data.plot <- melt(mean.indiv.gene, c(as.character(pars$meta_week),as.character(pars$meta_group)))

    #Convert ProbeID in illumina ID
    probes2ill <- as.list(illuminaHumanv4ARRAYADDRESS)
    ill2symb <- as.list(illuminaHumanv4SYMBOL)
    gene_name <- NULL
    for(i in 1:nrow(mean.indiv.data.plot)) {
      gene_name <-  c(gene_name,as.character(ill2symb[as.character(names(probes2ill[which(probes2ill %in% mean.indiv.data.plot[i,"variable"])]))]))
    }
    mean.indiv.data.plot$Gene <- gene_name

    #Normalize
    mean.indiv.data.plot$Norm <- rep(NA,nrow(mean.indiv.data.plot))
    for(gene in unique(mean.indiv.data.plot$Gene)) {
      if(!is.null(meta_group)) {
        for(group in unique(mean.indiv.data.plot[,as.character((pars$meta_group))])) {
          index <- which(mean.indiv.data.plot$Gene %in% gene & unique(mean.indiv.data.plot[,as.character((pars$meta_group))]) %in% group)
          mean.indiv.data.plot[index,]$Norm <- normal_distribution(mean.indiv.data.plot[index,]$value)
        }
      }
      else {
        index <- which(mean.indiv.data.plot$Gene %in% gene)
        mean.indiv.data.plot[index,]$Norm <- normal_distribution(mean.indiv.data.plot[index,]$value)
      }
    }

    mean.indiv.data.plot[,as.character(pars$meta_week)] <- as.numeric(as.character(mean.indiv.data.plot[,as.character(pars$meta_week)]))
    #Plot
    colnames(mean.indiv.data.plot)[which(colnames(mean.indiv.data.plot) %in% as.character(pars$meta_week))] <- "week"
    if(!is.null(meta_group)) colnames(mean.indiv.data.plot)[which(colnames(mean.indiv.data.plot) %in% as.character(pars$meta_group))] <- "group"
    p <- ggplot(data=mean.indiv.data.plot,aes(x=week, y = Norm, colour = Gene), na.rm = TRUE) +
      geom_point(aes(group=Gene), size=1) +
      geom_line(aes(group=Gene), linetype='dashed') +
      ylab(label = "Gene expression")
      if(!is.null(meta_group)) p <- p + facet_wrap(~ group, ncol = 2)
    p <- p + ggtitle("Dynamic of gene expression in each arm") +
      geom_vline(xintercept=unique(mean.indiv.data.plot$week),linetype=4, color="#A8A8A8") +
      scale_x_continuous(breaks = unique(mean.indiv.data.plot$week)) +
      theme_bw()

    print(p)

    return(mean.indiv.data.plot)

    ##Return plot and data frame
  }

  ##Option with sex and age and group
  else {
    ###TO DO
    #Return plot and data frame
    print("This part is in construction...")
  }


}
