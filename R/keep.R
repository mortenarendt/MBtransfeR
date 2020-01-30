#' empty the workspace keeping only the listed objects
#'
#'
#'@export
keep <- function(x = c()){
  lss <- ls()
  lss <- c(lss[lss!=x],'lss')
  rm(list = lss)
}