#' Calculate transfer statistics including permutation testing for two phyloseq objects
#' 
#'@param phy1 phyloseq object for the giving environment
#'@param phy2 phyloseq object for the recieving environment
#'@param identifier  in sample_data for mathcing dyads(character)
#'@import dplyr tidyverse broom
#'@export
randpermutationTransferStats <- function(phy1, phy2, identifier, nperm = 1000, verbose = F){
  # dig out the count_tables 
  TAXtb <- data.frame(tax_table(phy1))
  o1 <- phy1 %>% otu_table()
  o2 <- phy2 %>% otu_table()
  if (taxa_are_rows(phy1)){
    o1 <- o1 %>% t()
    o2 <- o2 %>% t()
  } 
  # find intersect of the two datasets
  ID1 <- sample_data(phy1)[[identifier]]
  ID2 <- sample_data(phy2)[[identifier]]
  
  # order according to ID1 and ID2
  id1 <- order(ID1)
  id2 <- order(ID2)
  ID1 <- ID1[id1]
  ID2 <- ID2[id2]
  o1 <- o1[id1,]
  o2 <- o2[id2,]
  
  # find intersect
  IDmatch <- intersect(ID1,ID2)
  
  # cut out samples
  o1 <- o1[ID1 %in% IDmatch,]
  o2 <- o2[ID2 %in% IDmatch,]
  ID1 <- ID1[ID1 %in% IDmatch]
  ID2 <- ID2[ID2 %in% IDmatch]
  
  deptho1 <- apply(o1,1,sum)
  deptho2 <- apply(o2,1,sum)
  totAbuM <- sum(deptho1)
  totAbuC <- sum(deptho2)
  
  # remove non-testabel OTU'
  s01 <- apply(o1==0,2,sum)
  s02 <- apply(o2==0,2,sum)
  nn <- dim(o2)[1]
  icvar <- s01>0 & s02>0 & s01<nn & s02<nn
  o1 <- o1[,icvar]
  o2 <- o2[,icvar]
  p <- sum(icvar)
  print(dim(o1))
  TAXtb <- TAXtb[icvar,]
  
  abuM <- apply(o1,2,sum)
  abuC <- apply(o2,2,sum)
  abuMrel <- o1 %>% sweep(1,deptho1,'/') %>% apply(2,mean)
  abuCrel <- o2 %>% sweep(1,deptho2,'/') %>% apply(2,mean)
  
  STAT <- get2by2table(o1,o2)
  STAT <- data.frame(STAT,abuM,abuMrel,abuC,abuCrel)
  STAT$otu <- rownames(STAT)
  
  # make permutation testing
  if (nperm>0){
    permSTAT <- array(dim = c(dim(STAT)[1],4,nperm))
    for (i in 1:nperm){
      stp <- get2by2table(o1[sample(dim(o1)[1]),],o2)
      permSTAT[,,i] <-  as.matrix(stp)
    }
    dimnames(permSTAT)[[1]] <- rownames(stp)
    dimnames(permSTAT)[[2]] <- colnames(stp)
    
    permSTATfisher <- array(dim = c(dim(STAT)[1],4,nperm))
    for (i in 1:nperm){
      if (verbose) {print(i)}
      aa <- permSTAT[,,i] %>% data.frame()
      for (j in 1:dim(permSTAT)[1]){
        stp2 <- getFisher(aa[j,],doglm = F)
        # permSTATfisher[j,,i] <- c(stp2$Fisher_estimate, stp2$Fisher_p.value)
        permSTATfisher[j,,i] <- c(stp2$Fisher_estimate,stp2$or_biascorr, stp2$Gtest_p.value, stp2$Fisher_p.value)
      }
    }
    dimnames(permSTATfisher)[[1]] <- rownames(stp)
    # dimnames(permSTATfisher)[[2]] <- c('Fisher_estimate', 'Fisher_p.value')
    dimnames(permSTATfisher)[[2]] <- c('Fisher_estimate','or_biascorr', 'Gtest_p.value', 'Fisher_p.value')
  } 
  else {
    permSTAT <- NULL
    permSTATfisher = NULL
  }
  
  STAT <- STAT %>%
    group_by(otu) %>%
    do(getFisher(x = .)) %>%
    ungroup %>%
    left_join(STAT,by = 'otu')
  
  STAT <- merge(STAT,TAXtb,by.x = 'otu',by.y = 'row.names')
  STAT$totAbuM <- totAbuM
  STAT$totAbuC <- totAbuC  
  
  STAT <- STAT %>%
    mutate(Fisher_estimatetr =  truncateZerosInf(Fisher_estimate),
           Glm_ortr =  truncateZerosInf(Glm_or))
  
  return(list(STAT,permSTAT, permSTATfisher))
} 