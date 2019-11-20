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
n_list <- c(10,50,100,500,1000,5000,10000)
pi <- 1
alpha <- 1
beta <- 1
alpha=0.05
precBits = 1024*16

critical_list <- foreach(i=1:(3*length(n_list)),.errorhandling='pass',.packages="jointTest")%dopar%{
  simRep=100000
  k=(i-1)%/%length(n_list)
  i=i-length(n_list)*k
  if(k==0){
    MPlusCriticalMC(alpha,n_list[i],rep=simRep)
  }else{
    if(k==1){
      MMinusCriticalMC(alpha,n_list[i],rep=simRep)
    }else{
      BJCriticalMC(alpha,n_list[i],rep=simRep)
    }
  }
}
critical_list=matrix(critical_list,ncol=3)
load("typeI_critical")

mp_critical_list <- critical_list[,1]
mm_critical_list <- critical_list[,2]
bj_critical_list <- critical_list[,3]

result=c()
for(k in seq_along(n_list)){
  n<- n_list[k]
  message(n)
  mp_critical <- mp_critical_list[k]
  mm_critical <- mm_critical_list[k]
  bj_critical <- bj_critical_list[k]
  mp_asym_critical <- -log(1-alpha)/2/log(n)/log(log(n))
  mm_asym_critical <- -log(1-alpha)/2/log(n)/log(log(n))
  bj_asym_critical <- -log(1-alpha)/4/log(n)/log(log(n))
  
  record <- matrix(0,nRep,6)
  x <- matrix(
    runif(n*nRep),nRep,n)
  
  mp=apply(x,1,MPlusStat)
  mm=apply(x,1,MMinusStat)
  bj=apply(x,1,BJStat)
  
  record[,1] <- mp < mp_critical
  record[,2] <- mm < mm_critical
  record[,3] <- bj < bj_critical
  
  record[,4] <- mp < mp_asym_critical
  record[,5] <- mm < mm_asym_critical
  record[,6] <- bj < bj_asym_critical
  
  result <- rbind(result, colMeans(record))
}
#save(result,file="typeI_simulate")
load("typeI_simulate")

x=seq_along(n_list)
plot(x,result[,1],type="b",ylim=c(0,0.2),xaxt="n",xlab="Sample Size",ylab="Type I Error",lty = 1, lwd = 1,pch=1)
lines(x,result[,2],type="b",pch=2)
lines(x,result[,3],type="b",pch=3)

lines(x,result[,4],type="b",lty = 2,pch=1)
lines(x,result[,5],type="b",lty = 2,pch=2)
lines(x,result[,6],type="b",lty = 2,pch=3)

xlab=n_list
axis(side=1, at=seq_along(xlab), labels = xlab)


legend(x=6,y=0.22,
       legend=c("exact M","exact M+","exact M-","asym M","asym M+","asym M-"),
       lty=c(1,1,1,2,2,2),pch=rep(rep(1:3)),cex = 0.7,bty = "n"
       
       )


#as a function of alpha
nRep <- 10000
n <- 1000
pi <- 1
alpha_list <- seq(0,1,0.01)
beta <- 1
alpha=0.05
precBits = 1024*16



simRep=100000
mp_critical_list<-MPlusCriticalMC(alpha_list,n,rep=simRep)
mm_critical_list<-MMinusCriticalMC(alpha_list,n,rep=simRep)
bj_critical_list<-BJCriticalMC(alpha_list,n,rep=simRep)

#save(mp_critical_list,mm_critical_list,bj_critical_list,file="typeI_critical2")
load("typeI_critical2")


x <- matrix(
  runif(n*nRep),nRep,n)
mp=apply(x,1,MPlusStat)
mm=apply(x,1,MMinusStat)
bj=apply(x,1,BJStat)

result=c()
record <- matrix(0,nRep,6)
for(k in 1:length(alpha_list)){
alpha<- alpha_list[k]
mp_critical <- mp_critical_list[k]
mm_critical <- mm_critical_list[k]
bj_critical <- bj_critical_list[k]
mp_asym_critical <- -log(1-alpha)/2/log(n)/log(log(n))
mm_asym_critical <- -log(1-alpha)/2/log(n)/log(log(n))
bj_asym_critical <- -log(1-alpha)/4/log(n)/log(log(n))




record[,1] <- mp < mp_critical
record[,2] <- mm < mm_critical
record[,3] <- bj < bj_critical

record[,4] <- mp < mp_asym_critical
record[,5] <- mm < mm_asym_critical
record[,6] <- bj < bj_asym_critical

result <- rbind(result, colMeans(record))
}

#save(result,file="typeI_simulate2")
load("typeI_simulate2")

k=seq(0,101,10)


x=seq_along(alpha_list)
plot(alpha_list[k],result[k,1],type="b",xlab=expression(alpha),ylab="Estimated Type I Error",lty = 1, lwd = 1,pch=1)
lines(alpha_list[k],result[k,2],type="b",pch=2)
lines(alpha_list[k],result[k,3],type="b",pch=3)

lines(alpha_list[k],result[k,4],type="b",lty = 2,pch=1)
lines(alpha_list[k],result[k,5],type="b",lty = 2,pch=2)
lines(alpha_list[k],result[k,6],type="b",lty = 2,pch=3)


legend(x=0.84,y=0.6,
       legend=c("exact M","exact M+","exact M-","asym M","asym M+","asym M-"),
       lty=c(1,1,1,2,2,2),pch=rep(rep(1:3)),cex = 0.7,bty = "n"
       
)
