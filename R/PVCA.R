#' Performs a Principal Variance Component Analysis (PVCA) and output the plot
#' 
#' @param expression numeric matrix of gene expressions with sample IDs as \code{colnames} and gene/probes in rows
#' @param sample_info a data.frame containing the design and technical informations on the samples.
#' @param Sample_ID the name of the column in \code{sample_info} which contain the sample IDs and which should match \code{colnames(expression)} in the same order.
#' @param effects the name of the columns from \code{sample_info} which should be investigated for potential Batch effect
#' @param pct_threshold a value between 0 and 1 (Default is \code{0.95}). The percentile value of the minimum amount of the variabilities that the selected principal components need to explain.
#' @param verbose a logical indicating if information about method's progress must be printed. Default is \code{FALSE}.
#' @param title the title of the barplot.
#' @param ... others parameters to be passed into the \code{barplot}.
#' 
#' @return a barplot.
#' 
#' @importFrom lme4 lmer VarCorr
#' @importFrom graphics axis barplot par text
#' @importFrom stats cor na.omit sigma
#' 
#' @author Boris Hejblum
#' 
#' @references Adapted from Pierre R. Bushel original programm 
#' (email: \code{Bushel@niehs.nih.gov})
#' 
#' @export

PVCA <- function(expression, sample_info, Sample_ID, effects,
                 pct_threshold=0.95, 
                 verbose=FALSE, title="PVCA", ...){
  
  pars <- as.list(match.call()[-1])
  
  ########## Load data ##########
  
  if(!is.data.frame(sample_info)){
    sample_info <- as.data.frame(sample_info)
  }
  
  if(any(!(colnames(expression) %in% sample_info[, as.character(pars$Sample_ID)])) | 
     any(!sample_info[, as.character(pars$Sample_ID)] %in% colnames(expression))){
    stop("colnames(expression) and sample_info$Sample_ID 
         do not contain the same elements")
  }
  
  if(any(as.character(colnames(expression))!= sample_info[, as.character(pars$Sample_ID)])){
    stop("colnames(expression) and sample_info$Sample_ID do not match exactly.\n 
         They are probably ordered differently...")
  }
  
  if(any(colSums(is.na(sample_info[, effects])) > 0)){
    stop("missing values in column ", names(which(colSums(is.na(sample_info[, effects]))>0)), " from `sample_info`")
  }
  
  theDataMatrix <- expression[,order(as.character(colnames(expression)))]
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  
  sample_info <- sample_info[order(sample_info[, as.character(pars$Sample_ID)]),]
  expDesignRowN <- nrow(sample_info)
  expDesignColN <- ncol(sample_info)
  myColNames <- names(sample_info)
  
  # Amount of variability desired to be explained by the principal components. Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0
  
  ########## Center the data (center rows) ##########
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
  theDataMatrixCentered_transposed <-  apply(X=theDataMatrix, MARGIN=1, FUN=scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered <-  t(theDataMatrixCentered_transposed)
  
  ########## Compute correlation matrix ##########
  
  theDataCor <- cor(theDataMatrixCentered)
  
  ########## Obtain eigenvalues ##########
  
  eigenData <- eigen(theDataCor)
  eigenValues <- eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix <- eigenData$vectors
  eigenValuesSum <- sum(eigenValues)
  percents_PCs <- eigenValues/eigenValuesSum 
  
  ########## Merge experimental file and eigenvectors for n components ##########
  
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1){
    
    my_sum_2  = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= pct_threshold){
      my_counter_2 = my_counter_2 + 1
    }
    if(verbose){message(paste(my_counter_2,"      ", 1 - my_sum_2))}
  }
  if (my_counter_2 < 3){
    pc_n  = 3
    
  }else {
    pc_n <- my_counter_2 
  }
  
  # pc_n is the number of principal components to model
  explained_var <- sum(percents_PCs[1:pc_n])
  
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
  mycounter <- 0
  for (i in 1:pc_n){
    for (j in 1:expDesignRowN){
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter,1] <- eigenVectorsMatrix[j,i]
    }
  }
  
  AAA <- sample_info[rep(1:expDesignRowN,pc_n),]
  
  Data <- cbind.data.frame(AAA,pc_data_matrix)
  
  par(mfrow=c(1,1))
  
  effects_n <- length(effects) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  formula <- paste("pc_data_matrix ~ (1|`", effects[1],"`)",sep='')
  for (i in 2:length(effects)){
    formula <- paste(formula," + (1|`", effects[i],"`)",sep='')
  }

  y <- 0
  # browser()
  for (i in 1:pc_n){
    y <- (((i-1)*expDesignRowN)+1)
    
    if(verbose){
      message("fit PC #", i)
    }
   
    Rm1ML <- (suppressMessages(lmer(
      formula = formula,
      data = Data[y:(((i-1)*expDesignRowN) + expDesignRowN), ], REML = TRUE,  
      verbose = FALSE, na.action = na.omit)))
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), resid = sigma(Rm1ML)^2)
  }
  
  effectsa <- c(as.character(names(VarCorr(Rm1ML))), "residual")
  
  ########## Standardize Variance ##########
  
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    mySum = sum(randomEffectsMatrix[i,])
    for (j in 1:effects_n){
      randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
    }
  }
  
  ########## Compute Weighted Proportions ##########
  
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    weight <- eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n){
      randomEffectsMatrixWtProp[i,j] <- randomEffectsMatrixStdze[i,j]*weight
    }
  }
  
  ########## Compute Weighted Ave Proportions ##########
  
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp,na.rm=TRUE)
  totalSum <- sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
  
  for (j in 1:effects_n){
    randomEffectsMatrixWtAveProp[j] <- randomEffectsSums[j]/totalSum 	
  }
  
  message(paste(pc_n,' axes explaining ', formatC(explained_var*100, digits = 2), '% of variance', sep=''))
  
  bp <- barplot(randomEffectsMatrixWtAveProp,  xlab = "", ylab = "Weighted average proportion variance", 
                ylim= c(0,1.1),col = c("blue"), las=2, main=title,
                ...)
  axis(1, at = bp, labels = FALSE, xlab = "Effects", cex.axis = 1, las = 2)
  text(bp, par("usr")[3]+1, adj = 1,labels = effectsa, srt = 90, xpd = TRUE, cex = 1.2)
  values <- randomEffectsMatrixWtAveProp
  new_values <- round(values, 2)
  text(bp,randomEffectsMatrixWtAveProp, labels = new_values, pos = 3, cex = 1)
}