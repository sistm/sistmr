#' Genes dynmaic over time
#'
#' @param data data
#' @param meta
#' @param vec_week
#' @param vec_group
#' @param vec_gene
#' @param meta_week
#' @param meta_group
#' @param meta_sample_ID
#' @param meta_participant
#' @param meta_sex
#' @param meta_age
#' @param indiv
#' @param group_facet
#' @param convert
#' @param legend
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
#'
dynamic_genes <- function(data, meta, vec_week, vec_group=NULL,
                  vec_gene,meta_week, meta_group=NULL, norm_group=NULL,
                  meta_sample_ID,meta_participant,meta_sex=NULL, meta_age=NULL,
                  indiv,group_facet=FALSE, convert=FALSE, legend=TRUE,
                  path_output=NULL, nameFile=NULL, title=NULL,norm="reduce_center",
                  indiv_col=TRUE, group_col=FALSE, norm_ind=FALSE,indiv_return_plot=1,
                  indiv_return_data=0) {

  pars <- as.list(match.call()[-1])
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


  #Mean individus
  if(indiv == 0) {

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


    if(convert==TRUE) {
    probes2ill <- as.list(illuminaHumanv4ARRAYADDRESS)
    ill2symb <- as.list(illuminaHumanv4SYMBOL)
    gene_name <- NULL
    for(i in 1:nrow(mean_indiv_data_plot)) {
      gene_name <-  c(gene_name,as.character(ill2symb[as.character(names(
                    probes2ill[which(probes2ill %in% mean_indiv_data_plot
                    [i,"variable"])]))]))
    }
    mean_indiv_data_plot <- mean_indiv_data_plot[,
                                                 !names(mean_indiv_data_plot) %in% "variable"]
    mean_indiv_data_plot[,"variable"] <- gene_name
    }
    #################################
    #Normalize

    if(norm == "reduce_center") {
      mean_indiv_data_plot$Norm <- rep(NA,nrow(mean_indiv_data_plot))
      for(gene in unique(mean_indiv_data_plot$variable)) {
        if(!is.null(pars$meta_group) & !is.null(pars$norm_group)) {
          for(groupi in unique(mean_indiv_data_plot
            [,as.character((pars$meta_group))])) {
            index <- which(mean_indiv_data_plot$variable %in% gene &
                     mean_indiv_data_plot[,as.character((pars$meta_group))]
                      %in% groupi)
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
    }

    if(norm=="around_zero") {
      mean_indiv_data_plot$Norm <- rep(NA,nrow(mean_indiv_data_plot))
      for(gene in unique(mean_indiv_data_plot$variable)) {
        if(!is.null(pars$meta_group)) {
          for(groupi in unique(mean_indiv_data_plot
                               [,as.character((pars$meta_group))])) {
            index <- which(mean_indiv_data_plot$variable %in% gene &
                             mean_indiv_data_plot[,as.character((pars$meta_group))]
                              %in% groupi)
            mean_indiv_data_plot[index,]$Norm <- normal_zero(
              mean_indiv_data_plot[index,]$value)
          }
        }
        else {
          index <- which(mean_indiv_data_plot$variable %in% gene)
          mean_indiv_data_plot[index,]$Norm <-
            normal_zero(mean_indiv_data_plot[index,]$value)
        }
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

    p <- p + ggtitle(title) +
         geom_vline(xintercept=unique(mean_indiv_data_plot$time),
         linetype=4, color="#A8A8A8") +
         scale_x_continuous(breaks =
         unique(mean_indiv_data_plot$time)) +
         theme_bw()

    if(legend == FALSE) p <- p + theme(legend.position="none")

    return(list(p,mean_indiv_data_plot))

    ##Return plot and data frame
  }

  ##Option with sex and age and group

  else {
    ###TO DO
    #Open pdf
    list_p <- list()
    list_data <- list()
    if(is.null(path_output)) {
      path <- getwd()
    } else {
      path <- path_output
    }
    #pdf(file=paste0(path, "/", nameFile, ".pdf"))
    #For each genes
    i <- 1
    for(gene in vec_gene) {

        data_gene <- cbind(as.character(merge_data[,as.character(pars$meta_group)]),
                           as.character(merge_data[,as.character(pars$meta_week)]),
                           merge_data[,as.character(pars$meta_sample_ID)],
                           as.character(merge_data[,as.character(pars$meta_participant)]),
                           round(merge_data[,which(colnames(merge_data) %in% gene)],2))

        colnames(data_gene) <- c("group","time","sample_ID","participant","gene"
                                 )
        data_gene <- as.data.frame(data_gene)
        data_gene$Norm <- NA


        ##NORMALISATION
        if(!is.null(pars$meta_group) && norm_ind == TRUE) {
          if(norm_group == TRUE) {
            for(gr in unique(data_gene$group)) {
              for(ind in unique(data_gene$participant)) {
                if(norm=="reduce_center") {
                  data_gene[which(data_gene$participant %in% ind & data_gene$group %in% gr),"Norm"] <-
                  normal_distribution(data_gene[which(data_gene$participant %in% ind & data_gene$group %in% gr),"gene"])
                }
                if(norm=="around_zero") {
                    data_gene[which(data_gene$participant %in% ind & data_gene$group %in% gr),"Norm"] <-
                    normal_zero(data_gene[which(data_gene$participant %in% ind & data_gene$group %in% gr),"gene"])
                }

              }
            }
          }
          if(norm_group == FALSE) {
            for(ind in unique(data_gene$participant)) {
              if(norm=="reduce_center") {
                data_gene[which(data_gene$participant %in% ind),"Norm"] <-
                normal_distribution(data_gene[which(data_gene$participant %in% ind),"gene"])
              }
              if(norm=="around_zero") {
                data_gene[which(data_gene$participant %in% ind),"Norm"] <-
                normal_zero(data_gene[which(data_gene$participant %in% ind),"gene"])
              }
            }
          }
        }

        if(!is.null(pars$meta_group) && norm_ind == FALSE) {
          if(norm_group == TRUE) {
            for(gr in unique(data_gene$group)) {
                if(norm=="reduce_center") {
                  data_gene[which(data_gene$group %in% gr),"Norm"] <-
                    normal_distribution(data_gene[which(data_gene$group %in% gr),"gene"])
                }
                if(norm=="around_zero") {
                  data_gene[which(data_gene$group %in% gr),"Norm"] <-
                    normal_zero(data_gene[which(data_gene$group %in% gr),"gene"])
                }

            }
          }
          if(norm_group == FALSE) {
            for(ind in unique(data_gene$participant)) {
              if(norm=="reduce_center") {
                data_gene[,"Norm"] <-
                  normal_distribution(data_gene[,"gene"])
              }
              if(norm=="around_zero") {
                data_gene[,"Norm"] <-
                  normal_zero(data_gene[,"gene"])
              }
            }
          }
        }

        if(is.null(pars$meta_group) && norm_ind == TRUE) {
          for(ind in unique(data_gene$participant)) {
            if(norm=="reduce_center") {
              data_gene[which(data_gene$participant %in% ind),"Norm"] <-
                normal_distribution(data_gene[which(data_gene$participant %in% ind),"gene"])
            }
            if(norm=="around_zero") {
              data_gene[which(data_gene$participant %in% ind),"Norm"] <-
                normal_zero(data_gene[which(data_gene$participant %in% ind),"gene"])
            }
          }
        }

        if(is.null(pars$meta_group) && norm_ind == FALSE) {
            if(norm=="reduce_center") {
              data_gene[,"Norm"] <-
                normal_distribution(data_gene[,"gene"])
            }
            if(norm=="around_zero") {
              data_gene[,"Norm"] <-
                normal_zero(data_gene[,"gene"])
            }
        }



      data_gene$time <- as.numeric(as.character(data_gene$time))
        ##♣Plot
        if(group_facet == TRUE || is.null(pars$meta_group)) {
          p <- ggplot(data=data_gene,
                      aes(x=time, y = Norm,
                          colour = participant),na.rm = TRUE)
        }
        if(group_facet == TRUE || is.null(pars$meta_group) && group_col == TRUE) {
          p <- ggplot(data=data_gene,
                    aes(x=time, y = Norm,
                        colour = group),na.rm = TRUE)
        }
        if(group_facet == FALSE && !is.null(pars$meta_group)) {
          p <- ggplot(data=data_gene,
                      aes(x=time, y = Norm, colour = participant,
                          group = interaction(group, participant)),
                      na.rm = TRUE)

        }

        p <- p + geom_point(size=1) +
          geom_line(aes(group=participant),linetype='dashed') +
          ylab(label = "Gene expression")



        if(!is.null(pars$meta_group)) {
          if(group_facet == TRUE) p <- p + facet_wrap(~ group, ncol = 2)
        }
        p <- p + ggtitle(paste0("Dynamic of gene ",gene," expression")) +
          geom_vline(xintercept=unique(data_gene$time),
                     linetype=4, color="#A8A8A8") +
          scale_x_continuous(breaks =
                               unique(data_gene$time)) +
          theme_bw()

        if(legend == FALSE) p <- p + theme(legend.position="none")

        list_p[[i]] <- p
        list_data[[i]] <- data_gene
        names(list_data)[i] <- gene
        i <- i + 1


     }

      if(indiv_return_plot == 1) return(list_p)
      if(indiv_return_data == 1) return(list_data)

  }
}



