normal_distribution <- function(vec) {
  vec <- as.numeric(vec)
  vec[is.na(vec)] <- mean(vec[!is.na(vec)])
  norm <- (vec - mean(vec))/sd(vec)

  return(norm)
}

normal_zero <- function(vec) {
  vec <- as.numeric(vec)
  norm <- (vec - vec[1])

  return(norm)
}