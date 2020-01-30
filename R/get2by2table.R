#' Calculates all 2x2 presens absense tables between two similar sized, and matched matrices
#'@param o1 matrix of counts (n x p) - n is samples and p is variables
#'@param o2 matrix of counts (n x p)
#'@import dplyr
#'@return a data.frame (p x 4) with the 2x2 table Ns (stats)
#'@export
get2by2table <- function(o1,o2){
  n11 <- t((o1>0)+0) %*% (o2>0) %>% diag
  n00 <- t((o1==0)+0) %*% (o2==0) %>% diag
  n10 <- t((o1>0)+0) %*% (o2==0) %>% diag
  n01 <- t((o1==0)+0) %*% (o2>0) %>% diag
  STAT <- data.frame(n00,n10,n01,n11)
  return(STAT)
}
