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




nRep <- 10000
n <- 1000
pi <- 0.8
alpha <- 1
alpha0 <- 0.05
beta <-2
precBits = 1024*16

loop_list <- seq(0,1,by=0.02)[-1]


message(n)
load("critical")
mp_critical<-as.numeric(critical_list[5,1])
HC_critical<-HCCriticalMC(alpha0,n,100000)


result = foreach(i = 1:length(loop_list),.combine = rbind,.packages = "jointTest")%dopar%{
  message(i)
  pi <- loop_list[i]
  
  x <- matrix(BUM(n*nRep,pi,alpha,beta),nRep,n)
  
  mp=apply(x,1,MPlusStat)
  hc=apply(x,1,HCStat)
  
  mean(mp<mp_critical)
  mean(hc>HC_critical)
  ks=apply(x,1,function(x)ks.test(x,punif)$p.value)
  mean(ks<0.05)
  bj_rj_rate <- mean(mp<mp_critical)
  hc_rj_rate <- mean(hc>HC_critical)
  c(bj_rj_rate,hc_rj_rate)
}

plot(seq_along(loop_list),result[,1],type="l",xaxt="n",xlab="pi",ylab="Probability of Rejection")
lines(seq_along(loop_list),result[,2],type= "l",lty=2)
xlab=round(seq(0.1,1,by=0.1),digits=2)
ind=which(round(loop_list,digits=2)%in%xlab)
axis(side=1, at=seq_along(loop_list)[ind], labels = xlab)
legend(x="topright",
       legend=c("BJ","HC"),bty = "n",
       lty=c(1,2),cex=1
)


a=seq(0.01,1,by=0.01)[10]
b=seq(0.1,1,by=0.1)[1]
a
b
a==b




legend(x="topright",
       legend=c("One-sided BJ","HC"),bty = "n",
       lty=c(1,2),cex=0.75
)

save(loop_list,result,file="sim2Result2")