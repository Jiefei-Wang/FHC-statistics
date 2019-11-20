
criticalMC<-function(func, alpha, n, rep){
  totalSample = 10000000
  
  rep_single = ceiling(totalSample/n)
  if(rep_single>rep)rep_single=rep
  
  outter = ceiling(n*rep/totalSample)
  
  
  stats=rep(0,rep_single*outter)
  for(i in 1:outter){
    x=matrix(runif(n*rep_single),n,rep_single)
    stats[(1+(i-1)*rep_single):(i*rep_single)]= apply(x,2,func)
  }
  return(quantile(stats,alpha))
}

#' @export
BJCriticalMC<-function(alpha,n,rep=10000){
  criticalMC(BJStat,alpha,n,rep)
}


#' @export
MPlusCriticalMC<-function(alpha,n,rep=10000){
  criticalMC(MPlusStat,alpha,n,rep)
}

#' @export
MMinusCriticalMC<-function(alpha,n,rep=10000){
  criticalMC(MMinusStat,alpha,n,rep)
}

HCCriticalMC<-function(alpha,n,rep=10000){
  criticalMC(HCStat,1-alpha,n,rep)
}