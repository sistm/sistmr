normal_distribution <- function(vec) {
  vec <- as.numeric(vec)
  vec[is.na(vec)] <- mean(vec[!is.na(vec)])
  norm <- (vec - mean(vec))/sd(vec)
  
  return(norm)
}