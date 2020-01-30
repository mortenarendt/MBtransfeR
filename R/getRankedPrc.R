#'
#'
#'@param X1
#'@param X2
#'@param identifier
#'@param groupid
#'@param metainfo
#'
#'@export
getRankedPrc <- function(X1,X2,identifier,groupid,metainfo){
  idvars <- c(identifier, metainfo,groupid)
  mX2 <- X2 %>% 
    gather(otu,Ccount,-idvars) %>%
    select(otu,Ccount,dyadnb)
  
  mX1 <- X1 %>% 
    gather(otu,Mcount,-idvars) %>%
    group_by(identifier) %>%
    arrange(desc(Mcount)) %>%
    mutate(rnk = 1:n()) %>%
    
    left_join(mX2,by = c('otu',identifier)) %>%
    filter(!is.na(Ccount)) %>%
    group_by(groupid,rnk) %>%
    summarise(n = n(), 
              nC = sum(Ccount>0),
              nC1 = sum(Ccount>1),
              nM = sum(Mcount>0)) %>%
    ungroup() %>%
    filter(nM==n) 
  return(mX1)
}
