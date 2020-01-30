#' Geometric mean
#'
#' Calculates the geometric mean of a vector (x)
#' @param x Numeric vector
#' @example 
#' x <- 1:100
#' gm_mean(x)
#' @return geometric mean (scalar) of input
#' @export
gm_mean <- function(x){
  exp(mean(log(x)))
}