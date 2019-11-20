library(foreach)
library(doParallel)
library(ddst)
cl=makeCluster(8)
registerDoParallel(cl)
# 
# z = runif(80)
# res=ddst.uniform.test(z, compute.p=TRUE)
# res$p.value


AF<-function(x,k=1.5){
  1-(1-x)^k
}

BF<-function(x,k=1.5){
  res=rep(0,length(x))
  ind=x<=0.5
  res[ind]=2^(k-1)*x[ind]^k
  res[!ind]=1-2^(k-1)*(1-x[!ind])^k
  res
}

CF<-function(x,k=1.5){
  res=rep(0,length(x))
  ind=x<=0.5
  res[ind]=0.5-2^(k-1)*(0.5-x[ind])^k
  res[!ind]=0.5+2^(k-1)*(x[!ind]-0.5)^k
  res
}

AF_inv<-function(x,k=1.5){
  1-(1-x)^(1/k)
}
BF_inv<-function(x,k=1.5){
  res=rep(0,length(x))
  ind=x<=0.5
  res[ind]=(2^(1-k)*x[ind])^(1/k)
  res[!ind]=1-(2^(1-k)*(1-x[!ind]))^(1/k)
  res
}
CF_inv<-function(x,k=1.5){
  res=rep(0,length(x))
  ind=x<=0.5
  res[ind]=0.5-(2^(1-k)*(0.5-x[ind]))^(1/k)
  res[!ind]=0.5+(2^(1-k)*(x[!ind]-0.5))^(1/k)
  res
}


sampleA<-function(n,k){
  x=runif(n)
  AF_inv(x,k)
}
sampleB<-function(n,k){
  x=runif(n)
  BF_inv(x,k)
}
sampleC<-function(n,k){
  x=runif(n)
  CF_inv(x,k)
}



KS<-function(x){
  sx=t(apply(x,1,sort))
  ksSeq = 1:ncol(x)/n
  ksPlus=abs(sweep(sx,2,ksSeq,"-"))
  ksMinus=abs(sweep(sx,2,ksSeq-1/n,"-"))
  ksPlus= apply(ksPlus,1,max)
  ksMinus= apply(ksMinus,1,max)
  ks=pmax(ksPlus,ksMinus)
}

KSCritical<-function(n, alpha0,nRep){
  x=matrix(runif(n*nRep),nRep,n)
  stat=KS(x)
  quantile(stat,1-alpha0)
}


BJCriticalList<- foreach(n = loop_list,.combine = c,.packages = "jointTest")%dopar%{
  BJCriticalMC(alpha0,n)
}

KSCriticalList<- foreach(n = loop_list,.combine = c,.packages = "jointTest")%dopar%{
  KSCritical(n,alpha0,10000)
}



nRep <- 10000
n <- 100
k=1.5
precBits = 1024*16
alpha0=0.05

loop_list <- c(10,20,40,50,70,100,150,200,400,500)
  # seq(1,1.5,by=0.05)




message(n)


result = foreach(i = 1:length(loop_list),.combine = rbind,.packages = c("jointTest","ddst"))%dopar%{
  message(i)
  n <- loop_list[i]
  
  bj_critical <- BJCriticalList[i]
  ks_critical<-KSCriticalList[i]
  
  x <- matrix(sampleB(n*nRep,k),nRep,n)
  
  bj=apply(x,1,BJStat)
  ks=KS(x)
  NST=apply(x,1,function(x)ddst.uniform.test(x, compute.p=TRUE)$p.value)
  
  
  bj_rj_rate <- mean(bj<bj_critical)
  ks_rj_rate <- mean(ks>ks_critical)
  NST_rj_rate <- mean(NST<0.05)
  c(bj_rj_rate,ks_rj_rate,NST_rj_rate)
}

plot(seq_along(loop_list),result[,1],type="l",ylim = c(0,1),xaxt="n",xlab="n",ylab="Probability of Rejection")
lines(seq_along(loop_list),result[,2],type= "l",lty=2)
lines(seq_along(loop_list),result[,3],type= "l",lty=3)
axis(side=1, at=seq_along(loop_list), labels = loop_list)

legend(x="topleft",
       legend=c("BJ","KS","Neyman"),bty = "n",
       lty=c(1,2,3),cex=0.75
)



save(result,file="sim2_simulate")








k=1.5
x=seq(0,1,by=0.001)
y=CF(x,k)
ind= which.max(abs(x-y))

x[ind]

plot(x,y,type="l",ylab="Prob",lty=1)
abline(0,1,lty=2)
abline(v=x[ind],lty=3)

axis(3,at=x[ind],labels= paste0("x=",x[ind]))

legend(x="bottomright",
       legend=c("CF(x)","uniform"),bty = "n",
       lty=c(1,2),cex=0.75
)



which.max(abs(diff(AF(x)-x)))

data=sampleC(10000,1.5)
hist(data)
plot(ecdf(data))
plot(x,BF(x),'l')

hist(AF_inv(x))




x=sampleA(10000,k)
