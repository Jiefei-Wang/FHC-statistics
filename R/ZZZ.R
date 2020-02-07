#' @import S4Vectors
#' @import Rmpfr
#' @importFrom gmp chooseZ
NULL

.onDetach<-function(libpath){
  print(libpath)
  cache$libpath <- libpath
  #save(list = "cache", file = "data/cachedData.RData")
}