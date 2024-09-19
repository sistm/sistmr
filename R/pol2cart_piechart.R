#' Convert polar coordinates to Cartesian coordinates
#'
#' @param theta 
#' @param r 
#' @param deg 
#' @param recycle 
#'
#' @return a \code{matrix} object
#' 
#' @keywords internal
#' 
#' @author Claire Bauduin (adapted to Binbin Xu codes) - Last modifications by Mélanie Huchon


pol2cart <- function(theta, r, deg = FALSE, recycle = FALSE) {
  if (deg)
    theta <- theta * pi / 180
  if (length(r) > 1 && length(r) != length(theta) && !recycle)
    stop(
      "'r' vector different length than theta, if recycling 'r' values is desired 'recycle' must be TRUE"
    )
  xx <- r * cos(theta)
  yy <- r * sin(theta)
  ## Account for machine error in trig functions
  xx[abs(xx) < 2e-15] <- 0
  yy[abs(yy) < 2e-15] <- 0
  out <- cbind(xx, yy)
  colnames(out) <- c("x", "y")
  return(out)
}