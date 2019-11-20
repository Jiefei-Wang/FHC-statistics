library(foreach)
library(doParallel)
cl=makeCluster(8)
registerDoParallel(cl)

BUM <- function(n ,pi ,alpha ,beta = 1){
  n1 <- rbinom(1, n, pi)
  x <- rep(0, n)
  ind <- sample.int(n,n1)
  x[ind] <- runif(n1)
  x[-ind] <- rbeta(n-n1, alpha,beta)
  x
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



nRep <- 2000
n <- 1000
pi <- 0.5
alpha <- 1
alpha0 <- 0.05
beta <-1
precBits = 1024*16

loop_list <- seq(0.5,1.5,by=0.05)


message(n)
bj_critical <- BJCriticalMC(alpha0,10000)
bj_critical<-as.numeric(critical_list[5,3])
ks_critical<-KSCritical(n,alpha0,100000)


result = foreach(i = 1:length(loop_list),.combine = rbind,.packages = "jointTest")%dopar%{
  message(i)
  beta <- loop_list[i]
  
  x <- matrix(BUM(n*nRep,pi,alpha,beta),nRep,n)
  
  bj=apply(x,1,BJStat)
  ks=KS(x)
  
  bj_rj_rate <- mean(bj<bj_critical)
  ks_rj_rate <- mean(ks>ks_critical)
  c(bj_rj_rate,ks_rj_rate)
}

plot(seq_along(loop_list),result[,1],type="l",ylim = c(0,1.1),xaxt="n",xlab="beta",ylab="Probability of Rejection")
lines(seq_along(loop_list),result[,2],type= "l",lty=2)
axis(side=1, at=seq_along(loop_list), labels = loop_list)

legend(x="topright",
       legend=c("BJ","KS"),bty = "n",
       lty=c(1,2),cex=0.75
)

save(loop_list,result,file="sim2Result2")

















res=ks.test(x,punif)
res$p.value


