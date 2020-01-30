#' Function to extract p-values from a permutation test dataset
#'
#'@param permSTAT a three dimensional array with 1. = OTUs, 2. = statistics, 3. = random permutation iterations
#'@param modelratio 1x1 scalar with the observed modelratio
#'@param trm 1x1 scalar, at which level to truncate the OR (OR>trm := trm, OR<1/trm := 1/trm)
#'@return data.frame with one row including pv (permutation p.value), SElgratio (standar error based on zero-centered normal assumptions and representability from the permutation null distribution, and some central limit theorem), permmedian (the median of the permuted null distribution), modelratio (the model ratio - as specified by input)
#'
#'@export
extractPV <- function(permSTAT,modelratio,trm=100){
  niter <- dim(permSTAT)[3]
  tb <- c()
  for (i in 1:niter){
    xx <- permSTAT[,,i] %>%
      data.frame() %>%
      mutate(Fisher_estimatetr = truncateZerosInf(Fisher_estimate,trm)) %>%
      # do(getWeightedRatio(x = .))
      do(ratio = getWeightedRatio(Fisher_estimatetr,Fisher_p.value))
    tb <- rbind(tb,xx)
  }
  pv <- sum(tb$ratio>modelratio) / niter
  
  # estimate null distribution for ratio
  lgratio <- log(tb$ratio)
  SElgratio <- sqrt(sum(lgratio^2)/length(lgratio))
  
  #print(c(i,median(tb$ratio),modelratio))
  df <- data.frame(pv, SElgratio, permmedian = median(tb$ratio), modelratio)
  
  return(df)
}
