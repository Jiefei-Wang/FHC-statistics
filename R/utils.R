#load(file="data/cachedData.RData")
cache <- new.env()
cache$nList <- new.env()
cache$nList[["500"]] <- 1000L
getCache <- function(signature){
  cache[[signature]]
}
setCache <- function(signature,result){
  cache[[signature]] <- result
}



getSuggestedPrecision <- function(n){
  if(!is.null(cache$nList[[as.character(n)]])){
    return(cache$nList[[as.character(n)]])
  }
  nList <- as.numeric(names(cache$nList))
  closedNIndex <- which.min(abs(nList-n)) 
  closedN <- nList[closedNIndex]
  dist <- closedN-n
  cache$nList[[as.character(nList)]]-dist
}

addSuggestedPrecision<- function(n,prec){
  if(!is.null(cache$nList[[as.character(n)]])){
    cache$nList[[as.character(n)]]<- max(cache$nList[[as.character(n)]],prec)
  }else{
    cache$nList[[as.character(n)]]<- prec
  }
}

saveCache <- function(){
  save(list = "cache", 
       file = paste0(cache$libpath,"/data/cachedData.RData"))
}

#' @export
print.jointTest <- function(x,...){
  print(as.numeric(x))
  invisible(x)
}


.jointTest<-function(statName,statValue,n,alpha0,index){
  if(!myMissing(index)){
    attris <-list(
      stat = statName,
      n=n,
      index= index
    )
  }else{
    attris <-list(
      stat = statName,
      n=n,
      alpha0= alpha0
    )
  }
  attris$class <- "jointTest"
  attributes(statValue)<-attris
  statValue
}


getArgs<-function(stat,n,alpha0,index){
  if(myMissing(n)){
    if(!is.null(attr(stat,"n"))) n <- attr(stat,"n")
    else stop("The sample size `n` is myMissing")
  }
  if(myMissing(alpha0)&&myMissing(index)){
    alpha0 <- attr(stat,"alpha0")
    index <- attr(stat,"index")
    if(is.null(alpha0)&&is.null(index)){
      index <- 1:n
      alpha0 <-1
    }else if(is.null(index)){
      nRegion <- max(floor(alpha0*n),1)
      index <- seq(1,nRegion)
    }
  }else if(myMissing(index)){
    nRegion <- max(floor(alpha0*n),1)
    index <- seq(1,nRegion)
  }
  list(n=n,index=index)
}

getPValueSigniture <-function(stat,n,alpha0,index,precBits,method,type){
  sigList <- list()
  if(myMissing(alpha0)){
    sigList$index <- index
  }else{
    sigList$alpha0 <- alpha0
  }
  sigList$n <- n
  sigList$stat <- unclass(stat)
  sigList$method <- method
  sigList$type <- type
  sigList$precBits <- precBits
  signature <- digest::digest(sigList)
}

getIndex<-function(index,indexL,indexU){
  if(myMissing(index)){
    if(myMissing(indexU)&&myMissing(indexL)){
      stop("indexU and indexL must be specified when index is myMissing")
    }else{
      if(myMissing(indexL)) indexL <- integer(0)
      if(myMissing(indexU)) indexU <- integer(0)
      
      index <- list(indexU = indexU,
                    indexL= indexL)
    }
  }else{
    if(!myMissing(indexU)||!myMissing(indexL)){
      stop("indexU and indexL should be empty when index is specified.")
    }
    index <- list(indexU = index,
                  indexL= index)
  }
  index
}


myMissing<-function(x){
  tryCatch({
    isMissing <- TRUE
    x
    isMissing <- FALSE
  },error =function(e) e,warning = function(w)w)
  isMissing
}