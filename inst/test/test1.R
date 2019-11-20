library(tictoc)
tic()
FHCPvalue(0.001,500,precBits=64)
toc()


stat=0.01
n=55
precBits=1024
nRegion=max(n,1)
l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
l[nRegion+seq_len(n-nRegion)]=l[nRegion]



record=c()
S_last_col1=getSCol(l,n,64)
S_last_col=getSCol(l,n,precBits)
for(j in (n-1):1){
  S_col1=getSCol(l,j,64)
  S_col=getSCol(l,j,precBits)
  tmp1=-S_col1*S_last_col1[j+1]
  tmp=-S_col*S_last_col[j+1]
  S_last_col1[1:j]=S_last_col1[1:j]+tmp1
  S_last_col[1:j]=S_last_col[1:j]+tmp
  theDif=max(as.numeric(S_last_col1-S_last_col))
  record=rbind(record,c(j,theDif))
}

k=n-1-(80:90)
plot(record[k,1],record[k,2],ylim=c(0,10))

abs(as.numeric(S_last_col1[1]))
abs(as.numeric(S_last_col[1]))


getSCol<-function(l,j,precBits){
  i_range=1:j
  tmp=mpfr(max(1-l[j],0),precBits=precBits)
  inner=tmp^(j-i_range+1)
  chooseNum=mpfr(chooseZ(j,j-i_range+1),precB=precBits)
  return(chooseNum*inner)
}


record=c()
S_last_col=as.numeric(getSCol(l,n,64))
for(j in (n-1):1){
  if(j==55) stop()
  S_col1=getSCol(l,j,64)
  S_col=getSCol(l,j,precBits)
  tmp1=-S_col1*S_last_col1[j+1]
  tmp=-S_col*S_last_col[j+1]
  S_last_col1[1:j]=S_last_col1[1:j]+tmp1
  S_last_col[1:j]=S_last_col[1:j]+tmp
  theDif=max(as.numeric(S_last_col1)-as.numeric(S_last_col))
  record=rbind(record,c(j,theDif))
}


allowedError=0.01/n
j=1:n
j*(1-l)^j

res=c()
for(j in 1:(n-1)){
  S_col=getSCol(l,j,precBits)
  res[j]=S_col[1]
}

