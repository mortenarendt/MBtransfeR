#' Permutation function to generate null distribution of distances bewteen randomly matched microbiome, conditioning on the environment
#'
#'@param D12 a 2n x 2n distance matric matrix
#'@param dyad 2n x 1 vector of dyad identifiers
#'
#'@export
permDistpair <- function(D12,dyad, nperm = 100){
  D <- getDistPairs(D12,dyad)
  
  n <-dim(D12)[1]/2
  # scrample second half of D12
  Dp <- matrix(NA, nrow = nperm,ncol = n)
  for (i in 1:nperm){
    idp <- sample(1:n)
    dyadp <- dyad
    dyad[1:n] <- dyad[idp]
    Dp[i,] <- getDistPairs(D12,dyadp)
  }
  mDp <- apply(Dp,1,'mean') 
  sDp <- mean(apply(Dp,1,'sd') )
  mD <- mean(D)
  p_perm <- sum(mDp<mD)/length(mDp)
  res <- data.frame(meandist = mD, sdmeandist = sd(D), meanrandomdist = mean(mDp), 
             sdmeanrandomdist = sDp, pv = p_perm)
  return(res)
  
}
