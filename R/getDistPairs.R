#' function for calculating distances between dyads
#'
#'@param D12 a 2n x 2n distance matric matrix
#'@param dyad 2n x 1 vector of dyad identifiers
#'@export
getDistPairs <- function(D12,dyad){
  D <- c()
  c <- 0
  for (i in unique(dyad)){
    c <- c+1
    d <-D12[dyad==i,dyad==i]
    D[c] <- d[1,2]
  }
  return(D)
}
