#' Calculates the individual areas of the vulcano plot based on the OR and pv 
#'
#'@param or Numeric vector (length k) of odd ratios
#'@param pv Numeric vector  (length k) of p-values
#'@example
#'
#'@return dataset with np - number of positive odds (or>1), nn - number of negative odds (or<=1), ratio - the relative area between positive and negative odds overall  
#'@export
getWeightedRatio <- function(or,pv,e = 0.001){
  pv <- -log10(pv)
  np <- sum(or>=1)
  nn <- sum(or<1)
  # mass
  area <- log(or)*pv
  ratio <- (sum(area[area>=0]) + e)/(-sum(area[area<0]) + e)
  R <- data.frame(np,nn , ratio)
  return(R)
}