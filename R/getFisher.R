#' Calculation of OR based on Fisher exact tests and GLM
#'
#'@param x dataset obj (1 x 4) with the 2x2 cross table of named as n00, n01, n10, n11
#'@import dplyr broom DescTools
#'@example 
#'x <- data.frame(n00 = 10,n10 = 3, n01 = 2, n11 = 7)
#'getFisher(x)
#'@return (transfer) odds based on different methods, returned as a data.frame()
#'@export
getFisher <- function(x,doglm = TRUE){
  x$or <- (x$n00*x$n11) / (x$n10*x$n01)
  side <- ifelse(x$or<1,'less','greater')
  nn <- x$n00 + x$n10 + x$n01 + x$n11
  tb <- x[,c('n00','n01','n10','n11')] %>% 
    unlist() %>%
    matrix(2) %>%
    fisher.test(alternative = side) %>%
    tidy
  colnames(tb) <- paste('Fisher',colnames(tb),sep = '_')
  
  # bias corrected
  or_biascorr <- (x$n00 + 0.5)*(x$n11+0.5) / (x$n01 + 0.5)*(x$n10+0.5)
  t_biascorr <- (nn*(abs(x$n00*x$n11 - x$n10*x$n01) - 0.5*nn)^2) / ((x$n00 + x$n01)*(x$n00 + x$n10)*(x$n10 + x$n11)*(x$n01 + x$n11))
  pv_biascorr <- 2*(1 - pnorm(t_biascorr))
  biascorr <- data.frame(or_biascorr,t_biascorr,pv_biascorr)
  
  # Gtest accoring to 
  gtest <- x[,c('n00','n01','n10','n11')] %>% 
    unlist() %>%
    matrix(2) %>%
    DescTools::GTest() %>% 
    tidy
  colnames(gtest) <- paste('Gtest',colnames(gtest),sep = '_')
  tb <- cbind(tb,biascorr,gtest)
  
  # glm
  if (doglm){
    mom <- c(rep(0,x$n00), rep(0,x$n01), rep(1,x$n10),rep(1,x$n11))
    child <- c(rep(0,x$n00), rep(1,x$n01), rep(0,x$n10),rep(1,x$n11))
    mdl <- glm(child~mom,family = binomial('logit')) %>% 
      tidy %>%
      filter(term=='mom') %>%
      mutate(or = exp(estimate)) %>%
      select(-term)
    colnames(mdl) <- paste('Glm',colnames(mdl),sep = '_')
    tb <- cbind(tb,mdl)
  } 
  return(tb)
}
