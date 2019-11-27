#' @export
HCStat = function(p, alpha0=1){
  n= length(p)
  nRegion=max(floor(alpha0*n),1)
  sp = sort(p)
  sp[sp==0]=min(0.000001,sp[sp!=0])
  HC = sqrt(n)*(seq(1,n)/n-sp)/sqrt(sp*(1-sp))
  out = max(HC[seq_len(nRegion)])
  return(out)
}
#' @export
MPlusStat<-function(p, alpha0=1){
  n=length(p)
  nRegion=max(floor(alpha0*n),1)
  sp=sort(p)
  sp[sp==0]=min(0.000001,sp[sp!=0])
  FHCElement=sapply(seq_along(p),function(x)pbeta(sp[x],x,n-x+1))
  min(FHCElement[seq_len(nRegion)])
}
#' @export
MMinusStat<-function(p, alpha0=1){
  n=length(p)
  nRegion=max(floor(alpha0*n),1)
  sp=sort(p)
  sp[sp==0]=min(0.000001,sp[sp!=0])
  FHCElement=1-sapply(seq_along(p),function(x)pbeta(sp[x],x,n-x+1))
  min(FHCElement[seq_len(nRegion)])
}
#' @export
BJStat<-function(p, alpha0=1){
  n=length(p)
  nRegion=max(floor(alpha0*n),1)
  sp=sort(p)
  sp[sp==0]=min(0.000001,sp[sp!=0])
  M_plus=sapply(seq_along(p),function(x)pbeta(sp[x],x,n-x+1))
  M_minus=1-M_plus
  min(c(M_plus[seq_len(nRegion)],M_minus[seq_len(nRegion)]))
}
KSStat<-function(p,alpha0=1){
  n=length(p)
  sp=sort(p)
  nRegion=max(floor(alpha0*n),1)
  
  ksSeq = 1:n/n
  ksPlus=sp-(ksSeq-1/n)
  ksMinus= ksSeq-sp
  max(c(ksPlus[seq_len(nRegion)],ksMinus[seq_len(nRegion)]))
}

getC = function(x, a){
  out=(x+(a^2-a*(a^2+4*(1-x)*x)^0.5)/2)/(1+a^2);
  return(out)
}

getSCol<-function(l,m,j,precBits){
  i_range=1:j
  tmp=mpfr(pmax(m[i_range]-l[j],0),precBits=precBits)
  inner=tmp^(j-i_range+1)
  chooseNum=mpfr(chooseZ(j,j-i_range+1),precB=precBits)
  return(chooseNum*inner)
}

orderedProb<-function(l,m,precBits=1024,progress=FALSE){
  n=length(l)
  #for(i in 1:(n-1)){
  #  l[l[i]>l[(i+1):n]]=l[i]
  #}
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
  mymethod= 1-abs(as.numeric(S_last_col[1]))
  return(mymethod)
}


#' @export
HCPvalue<-function(stat,n,alpha0=1,precBits=1024,progress=FALSE){
  nRegion=max(floor(alpha0*n),1)
  c_const=stat/sqrt(n)
  c_nonConst=1:n/n
  l=getC(c_nonConst,c_const)
  l[nRegion+seq_len(n-nRegion)]=l[nRegion]
  res=orderedProb(l,rep(1,n),precBits,progress)
  if(res>1||res<0) 
    stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  res
}

#' @export
MPlusPvalue<-function(stat,n,alpha0=1,precBits=1024,progress=FALSE,autoPrecision = TRUE){
  nRegion=max(floor(alpha0*n),1)
  l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
  l[nRegion+seq_len(n-nRegion)]=l[nRegion]
  res=orderedProb(l,rep(1,n),precBits,progress)
  if(res>1||res<0){
    if(autoPrecision)
      MPlusPvalue(stat,n,alpha0,precBits*2,progress,autoPrecision)
    else
      stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  }else{
    res
  }
}
#' @export
MMinusPvalue<-function(stat,n,alpha0=1,precBits=1024,progress=FALSE,autoPrecision = TRUE){
  nRegion=max(floor(alpha0*n),1)
  m=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
  m[nRegion+seq_len(n-nRegion)]=1
  res=orderedProb(rep(0,n),m,precBits,progress)
  if(res>1||res<0){
    if(autoPrecision)
      MMinusPvalue(stat,n,alpha0,precBits*2,progress,autoPrecision)
    else
      stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  }else{
    res
  }
}
#' @export
BJPvalue<-function(stat,n,alpha0=1,precBits=1024,progress=FALSE,autoPrecision = TRUE){
  nRegion=max(floor(alpha0*n),1)
  l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
  m=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
  
  l[nRegion+seq_len(n-nRegion)]=l[nRegion]
  m[nRegion+seq_len(n-nRegion)]=1
  res=orderedProb(l,m,precBits,progress)
  if(res>1||res<0){
    if(autoPrecision)
      BJPvalue(stat,n,alpha0,precBits*2,progress,autoPrecision)
    else
      stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  }else{
    res
  }
}
#' @export
KSPValue<-function(stat,n,alpha0=1,precBits=1024,progress=FALSE,autoPrecision = TRUE){
  nRegion=max(floor(alpha0*n),1)
  l= -stat + 1:n/n
  m= 1:n/n-1/n + stat
  
  l[l<0]=0
  m[m>1]=1
  
  l[nRegion+seq_len(n-nRegion)]=l[nRegion]
  m[nRegion+seq_len(n-nRegion)]=1
  res=orderedProb(l,m,precBits,progress)
  if(res>1||res<0){
    if(autoPrecision)
      KSPValue(stat,n,alpha0,precBits*2,progress,autoPrecision)
    else
      stop("An numeric overflow has been found, please consider to increase the precision. result: ", res)
  }else{
    res
  }
}




#' @export
HCCritical<-function(alpha,n,alpha0=1,precBits=1024,searchRange=c(0,100),...){
  rootFunc=function(stat,...) HCPvalue(stat,n,alpha0=alpha0,precBits=precBits,...)-alpha
  res=uniroot(rootFunc,searchRange,extendInt="yes",...)
  res$root
}

#' @export
MPlusCritical<-function(alpha,n,alpha0=1,precBits=1024,...){
  rootFunc=function(stat,...) MPlusPvalue(stat,n,alpha0=alpha0,precBits=precBits,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}
#' @export
MMinusCritical<-function(alpha,n,alpha0=1,precBits=1024,...){
  rootFunc=function(stat,...) MMinusPvalue(stat,n,alpha0=alpha0,precBits=precBits,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}
#' @export
BJCritical<-function(alpha,n,alpha0=1,precBits=1024,...){
  rootFunc=function(stat,...) BJPvalue(stat,n,alpha0=alpha0,precBits=precBits,...)-alpha
  res=uniroot(rootFunc,c(0,1),...)
  res$root
}


computeSMatrix <- function(l,m){
  n=length(l)
  s=matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if(j-i+1>=0){
        s[i,j]=choose(j,j-i+1)*
          (m[i]-l[j])^(j-i+1)
      }
    }
  }
  s
}

approxFactorial<-function(n){
  n*log(n)-n
}


computeLSMatrix <- function(l,m){
  n=length(l)
  S=matrix(0,n,n)
  
  for(i in 1:n){
    #message(i)
    j=i:n
    tmp=m[i]-l[j]
    tmp[tmp<=0]=1
    
    S[i,j]=lchoose(j,j-i+1)+(j-i+1)*log(tmp)*(tmp!=0)
  }
  S
}




