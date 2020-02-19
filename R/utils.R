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
  cache$nList[[as.character(closedN)]]-dist
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
  if(!is.null(index)){
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


getArgs<-function(stat,n,alpha0,index,indexL=NULL,indexU=NULL){
  if(is.null(n)){
    if(class(stat)!="jointTest"){
      stop("The sample size `n` is missing")
    }
    n <- attr(stat,"n")
  }
  level1 <- !is.null(alpha0)||!is.null(index)
  level2 <- !is.null(indexL)||!is.null(indexU)
  if(level1||
     level2){
    if(level1&&level2){
      stop("Either alpha0 and index or indexL and indexU must be NULL")
    }else{
    if(level1){
      if(!is.null(alpha0)&&!is.null(index)){
        stop("Either alpha0 or index should be NULL")
      }else{
        if(!is.null(alpha0)){
          nRegion <- max(floor(alpha0*n),1)
          index <- seq(1,nRegion)
        }
      }
    }else{
      index <- list(indexL= indexL,indexU=indexU)
    }
    }
  }else{
    if(class(stat)!="jointTest"){
      alpha0 <- 1
      index <- seq_len(n)
    }else{
      index <- attr(stat,"index")
      alpha0 <- attr(stat,"alpha0")
      if(is.null(index)){
        nRegion <- max(floor(alpha0*n),1)
        index <- seq(1,nRegion)
      }
    }
  }
  
  list(n=n,alpha0 = alpha0, index=index)
}

getPValueSigniture <-function(stat,args,precBits,method,type){
  if(!is.null(args$alpha0)){
    args$index <- NULL
    args$indexL <- NULL
    args$indexU <- NULL
  }
  args$stat <- as.numeric(stat)
  args$method <- method
  args$type <- type
  args$precBits <- precBits
  signature <- digest::digest(args)
}

getIndex<-function(n,alpha0, index,indexL,indexU){
  if(is.null(index)){
    if(is.null(indexU)&&is.null(indexL)){
      nRegion <- max(floor(alpha0*n),1)
      tmp <- seq(1,nRegion)
      index <- list(indexU = tmp,
                    indexL= tmp)
    }else{
      index <- list(indexU = indexU,
                    indexL= indexL)
    }
  }else{
    if(!is.null(indexU)||!is.null(indexL)){
      stop("indexU and indexL should be empty when index is specified.")
    }
    index <- list(indexL = index,
                  indexU= index)
  }
  index
}
