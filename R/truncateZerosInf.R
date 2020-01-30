#' Truncate large or small (abs(log(or)) >>> 0) odds ratios to trc or 1/trc
#'@param or Numeric vector
#'@example
#'or <- 10^rnorm(10)
#'truncateZerosInf(or,trc = 5)
#'@return a vector of truncated odds ratios
#'@export 
truncateZerosInf <- function(or,trc = 100){
  trcm <- 10^-log10(trc)
  ornew <- or
  ornew[(is.infinite(or) & or>1) | or>trc ] <- trc
  ornew[(or==0 & or<1) | or<trcm] <- trcm
  return(ornew)
}
