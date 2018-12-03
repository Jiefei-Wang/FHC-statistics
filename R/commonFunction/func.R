library(mvtnorm)


#Compute the p-values for the t-test
computeP=function(x,y,p=NA){
  if(!is.na(p[1])){
    return(p)
  }
  size=dim(x)
  pvalue=matrix(NA,1,size[2])
  for(i in 1:size[2]){
    pvalue[i]=t.test(x[,i],y[,i])$p.value
  }
  return(pvalue)
}

#Compute the p-values for the t-test
#More efficient than the above one
computeP2=function(x,y,p=NA){
  if(!is.na(p[1])){
    return(p)
  }
  x.size=dim(x)
  y.size=dim(y)
  x.m=colMeans(x)
  y.m=colMeans(y)
  df=x.size[1]+y.size[1]-2
  s=(colSums(sweep(x, 2, x.m, `-`)^2)+colSums(sweep(y, 2, y.m, `-`)^2))/(df)
  
  t=abs((x.m-y.m)/sqrt(s*(1/x.size[1]+1/y.size[1])))
  pvalue=matrix(pt(t[],df,lower.tail = F)*2,1,x.size[2])
  
  return(pvalue)
}





#compute the covariance matrix
#Parms:
#n: number of block
#type: The type of covariance matrix, can be auto regression(auto), compound symmetric(cs), independent(ind)
#and blocked compound symmetric(block_cs).
#parm: The parameter required by the covariance matrix
getCovMat<-function(n,type,parm){
  
  covMat=matrix(0,n,n)
  if(type=="auto"){
    ro=parm
    for(i in 1:n)
      for(j in 1:n){
        covMat[i,j]=ro^(abs(i-j))
      }
  }
  if(type=="cs"){
    sigma1=parm[1]
    sigma2=parm[2]
    for(i in 1:n)
      for(j in 1:n){
        if(i==j){
          covMat[i,j]=sigma1+sigma2
        }else{
          covMat[i,j]=sigma2
        }
      }
  }
  if(type=="ind"){
    covMat=diag(rep(1,n))
  }
  
  if(type=="block_cs"){
    #sigma1,sigma2,block size
    sigma1=parm[1]
    sigma2=parm[2]
    block=rmultinom(1, n-parm[3], rep(1/parm[3],parm[3]))+1
    rowind=1
    for(i in 1:parm[3]){
      covMat[rowind:(rowind+block[i]-1),rowind:(rowind+block[i]-1)]=getCovMat(block[i],"cs",c(sigma1,sigma2))
      rowind=rowind+block[i]
    }
  }
  
  return(covMat)
}
