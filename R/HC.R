## Indivial level function
HCLevel<-function(index,n,x,sp){
  sqrt(n)*(seq(1,n)[index]/n-sp[index])/sqrt(sp[index]*(1-sp[index]))
}
BJPlusLevel<-function(index,n,x,sp){
  if(length(index)==0)return(numeric(0))
  sapply(seq_along(x)[index],function(x)pbeta(sp[x],x,n-x+1))
}
BJMinusLevel<-function(index,n,x,sp){
  if(length(index)==0)return(numeric(0))
  1-BJPlusLevel(index,n,x,sp)
}
BJLevel<-function(index,n,x,sp){
  indexU <- index$indexU
  indexL <- index$indexL
  c(BJPlusLevel(indexU,n,x,sp),BJMinusLevel(indexL,n,x,sp))
}
KSPlusLevel<-function(index,n,x,sp){
  seq(1,n)[index]/n-sp[index]
}
KSMinusLevel<-function(index,n,x,sp){
  sp[index] - (seq(1,n)[index]-1)/n
}
KSLevel<-function(index,n,x,sp){
  c(KSPlusLevel(index,n,x,sp),KSMinusLevel(index,n,x,sp))
}

## These functions return a set of level stat
partialLevelStat <- function(statFunc,x,alpha0,index){
  if(is.null(index)){
    nRegion <-max(floor(alpha0*length(x)),1)
    index <- seq(1,nRegion)
  }
  n <- length(x)
  sp <- sort(x)
  sp[sp==0] <- min(10^-6,sp[sp!=0])
  sp[sp==1] <- max(1-10^-6,sp[sp!=1])
  statFunc(index=index ,n = n,x = x,sp =sp)
}

## Statistics
#' @export
HCStat<-function(x,alpha0 = 1, index=NULL){
  stat <- max(partialLevelStat(statFunc = HCLevel,
                       x = x,
                       alpha0 = alpha0,
                       index= index
                       ))
  .jointTest("HC",stat,length(x),alpha0,index)
}

#' @export
BJStat<-function(x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL){
  index <- getIndex(length(x),alpha0, index,indexL,indexU)
  stat <- min(partialLevelStat(statFunc = BJLevel,
                       x = x,
                       alpha0 = alpha0,
                       index= index
  ))
  .jointTest("BJ",stat,length(x),alpha0,index)
}

#' @export
BJPlusStat<-function(x,alpha0 = 1, index=NULL){
  BJStat(x=x,alpha0=alpha0,indexU=index)
}
#' @export
BJMinusStat<-function(x,alpha0 = 1, index=NULL){
  BJStat(x=x,alpha0=alpha0,indexL=index)
}

#' @export
KSStat<-function(x,alpha0=1,index=NULL){
  stat <- max(partialLevelStat(statFunc = KSLevel,
                       x = x,
                       alpha0 = alpha0,
                       index= index
  ))
  .jointTest("KS",stat,length(x),alpha0,index)
}
#' @export
KSPlusStat <- function(x,alpha0=1,index=NULL){
  stat <- max(partialLevelStat(statFunc = KSPlusLevel,
                       x = x,
                       alpha0 = alpha0,
                       index= index
  ))
  .jointTest("KS+",stat,length(x),alpha0,index)
}

#' @export
KSMinusStat <- function(x,alpha0=1,index=NULL){
  stat <- max(partialLevelStat(statFunc = KSMinusLevel,
                       x = x,
                       alpha0 = alpha0,
                       index= index
  ))
  .jointTest("KS-",stat,length(x),alpha0,index)
}



getSCol<-function(l,m,j,precBits){
  i_range=1:j
  tmp=mpfr(pmax(m[i_range]-l[j],0),precBits=precBits)
  inner=tmp^(j-i_range+1)
  chooseNum=mpfr(chooseZ(j,j-i_range+1),precB=precBits)
  return(chooseNum*inner)
}

orderedProb<-function(l,m,indexL,indexU,precBits,progress,autoPrecision){
  if(length(indexL)!=0){
    l[-indexL] <- 0
  }else{
    l=rep(0,length(l))
  }
  if(length(indexU)!=0){
    m[-indexU] <- 1
  }else{
    m=rep(1,length(m))
  }
  n <- length(l)
  for(i in seq_len(n-1)){
    if(l[i]>l[i+1]) l[i+1] <- l[i]
    j <- n-i
    if(m[j] > m[j+1]) m[j] <- m[j+1]
  }
  if(any(l>m))return(0)
  if(autoPrecision){
    precBits<- getSuggestedPrecision(n)
  }
  
  S_last_col=getSCol(l,m,n,precBits)
  if(length(S_last_col)!=1){
    if(progress) pb <- txtProgressBar(min=0,max=n)
    for(j in (n-1):1){
      S_col=getSCol(l,m,j,precBits)
      tmp=-S_col*S_last_col[j+1]
      S_last_col[1:j]=S_last_col[1:j]+tmp
      if(progress) setTxtProgressBar(pb, n-j)
    }
    if(progress) close(pb)
  }
  mymethod= abs(as.numeric(S_last_col[1]))
  if(mymethod>1||mymethod<0){
    if(autoPrecision)
      return(
        orderedProb(l,m,indexL,indexU,ceiling(precBits*1.2),progress,autoPrecision))
    else
      stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  }else{
    if(autoPrecision){
      addSuggestedPrecision(n,precBits)
    }
    return(mymethod)
  }
}

getC = function(x, a){
  out=(x+(a^2-a*(a^2+4*(1-x)*x)^0.5)/2)/(1+a^2);
  return(out)
}
#' @export
HCPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,
                   precBits=1024,progress=FALSE,autoPrecision = TRUE){
  args <- getArgs(stat,n,alpha0,index)
  signature <- getPValueSigniture(stat,args,precBits,"HC","pvalue")
  cachedResult <- getCache(signature)
  if(!is.null(cachedResult)) return(cachedResult)
  
  
  n <- args$n
  index <- args$index
  c_const=stat/sqrt(n)
  c_nonConst=1:n/n
  l=getC(c_nonConst,c_const)
  m=rep(1,n)
  res=1-orderedProb(l,m,index,index,precBits,progress,autoPrecision)
  setCache(signature,res)
  res
}

#' @export
BJPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,
                       precBits=1024,progress=FALSE,autoPrecision = TRUE){
  BJPvalue(stat = stat, n=n, alpha0=alpha0,
           indexL= index, precBits = precBits,
           progress = progress, autoPrecision=autoPrecision)
}
#' @export
BJMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,
                        precBits=1024,progress=FALSE,autoPrecision = TRUE){
  BJPvalue(stat = stat, n=n, alpha0=alpha0,
           indexU= index, precBits = precBits,
           progress = progress, autoPrecision=autoPrecision)
}
#' @export
BJPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,indexL=NULL,indexU=NULL,
                   precBits=1024,progress=FALSE,autoPrecision = TRUE){
  args <- getArgs(stat,n,alpha0,index,indexL,indexU)
  signature <- getPValueSigniture(stat,args,precBits,"BJ","pvalue")
  cachedResult <- getCache(signature)
  if(!is.null(cachedResult)) return(cachedResult)
  
  n <- args$n
  index <- args$index
  if(!is.list(index)){
    tmp <- list()
    tmp$indexL <- index
    tmp$indexU <- index
    index <- tmp
  }
  
  l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
  m=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
  
  res=1-orderedProb(l,m,index$indexL,index$indexU,precBits,progress,autoPrecision)
  setCache(signature,res)
  res
}
#' @export
KSPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,
                   precBits=1024,progress=FALSE,autoPrecision = TRUE){
  args <- getArgs(stat,n,alpha0,index)
  
  signature <- getPValueSigniture(stat,args,precBits,"KS","pvalue")
  cachedResult <- getCache(signature)
  if(!is.null(cachedResult)) return(cachedResult)
  
  n <- args$n
  index <- args$index
  
  l= 1:n/n - stat
  m= stat + 1:n/n-1/n
  
  l[l<0]=0
  m[m>1]=1
  res=1-orderedProb(l,m,index,index,precBits,progress,autoPrecision)
  setCache(signature,res)
  res
}




#' @export
HCCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL,precBits=1024,
                     autoPrecision = FALSE,searchRange=c(0,100),...){
  rootFunc=function(stat,...) 
    HCPvalue(stat,n,alpha0=alpha0,index=index,
             precBits=precBits,autoPrecision=autoPrecision,...)-alpha
  res=uniroot(rootFunc,searchRange,extendInt="yes",...)
  res$root
}

#' @export
BJPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL,precBits=1024,
                         autoPrecision = FALSE,...){
  rootFunc=function(stat,...) 
    BJPlusPvalue(stat,n,alpha0=alpha0,index=index,
                 precBits=precBits,autoPrecision=autoPrecision,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}
#' @export
BJMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL,precBits=1024,
                          autoPrecision = FALSE,...){
  rootFunc=function(stat,...) 
    BJMinusPvalue(stat,n,alpha0=alpha0,index=index,
                  precBits=precBits,autoPrecision=autoPrecision,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}
#' @export
BJCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL,
                     indexL=NULL,indexU=NULL,precBits=1024,
                     autoPrecision = FALSE,...){
  rootFunc=function(stat,...) 
    BJPvalue(stat,n,alpha0=alpha0,index=index,indexL = indexL,indexU=indexU,
             precBits=precBits,autoPrecision=autoPrecision,...)-alpha
  res=uniroot(rootFunc,c(0,0.5),...)
  res$root
}

#' @export
KSCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL,precBits=1024,
                     autoPrecision = FALSE,...){
  rootFunc=function(stat,...) 
    KSPvalue(stat,n,alpha0=alpha0,index=index,
             precBits=precBits,autoPrecision=autoPrecision,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}



