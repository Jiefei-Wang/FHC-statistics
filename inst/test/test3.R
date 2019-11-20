stat=0.01
n=100
precBits=1024
nRegion=max(n,1)
l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
l[nRegion+seq_len(n-nRegion)]=l[nRegion]

j=5

getLogSCol<-function(l,j){
  i_range=1:j
  tmp=logNumVector(max(1-l[j],0))
  inner=tmp^(j-i_range+1)
  chooseNum=logNumVectorFromLogValue(lchoose(j,j-i_range+1))
  return(chooseNum*inner)
}

as.numeric(getSCol(l,j,precBits))
a= getLogSCol(l,j)


S_last_col_log=getLogSCol(l,n)
for(j in (n-1):1){
  S_col=getLogSCol(l,j)
  tmp=-S_col*S_last_col_log[j+1]
  S_last_col_log[1:j]=S_last_col_log[1:j]+tmp
}
S_last_col_log[1]



a=logNumVector(1:3)
