
#Combine the result from the parallel computing
resultCmb<-function(result1, ...){
  result=Map(rbind, result1, ...)
  return(result)
}
#Compute the asymptotic critical value for the HC statistic
computeHCASYCV = function(N,alpha){
  bn = sqrt(2*log(log(N)))
  cn = 2*log(log(N))+.5*(log(log(log(N))) - log(4*pi))
  out = 1/bn*(cn-log(log(1/(1-alpha))))
  return(out)
}


#Compute HC_i
computeHClist = function(case,control,p=NA){
  if(is.na(p[1]))
    p=computeP(control,case)
  sp = sort(p)
  N= length(sp)
  HC = sqrt(N)*(seq(1,N)/N-sp)/sqrt(sp*(1-sp))
  #HC=HC^2
  return(HC)
}
#The Hc statistic, maximum of the HC_i
computeHC = function(case,control,p=NA){
  HC=computeHClist(case,control,p)
  out = max(HC)
  return(out)
}

#Compute the FHC statistic
computeFHC_ind<-function(case,control,p=NA){
  if(is.na(p[1])){
    p=computeP(control,case)
  }
  n=length(p)
  PHC=sort(p)
  out=rep(NA,n)
  for(i in 1:n){
    out[i]=pbeta(PHC[i],i,n+1-i)
  }
  
  return(min(out))
}


#Compute the FHC statistic
#This is a function specifically design for the permutation.
computeFHC_permute<-function(case,control,p=NA,parms){
  if(length(parms)==2){
    cdf.list=parms[[1]]
    index=parms[[2]]
  }else{
    cdf.list=parms
    index=rep(1,ncol(case))
  }
  if(is.na(p[1])){
    p=computeP(case,control)
  }
  sp=sort(p)
  ind=which(index!=0)
  out=rep(NA,length(p))
  for(i in 1:length(ind)){
    j=ind[i]
    out[i]=1-cdf.list[[j]](sp[j])
  }
  return(max(out,na.rm = T))
}


#Compute the KS statistic
computeKS<-function(case,control,p=NA){
  if(is.na(p[1])){
    p=computeP(case,control)
  }
  ks=ks.test(p,"punif")
  return(1-ks$p.value)
}


