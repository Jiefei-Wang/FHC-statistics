B<-function(n, i, x){
  pbeta(x,i,n-i+1)
}

G<-function(n,i,x){
  1-B(n,i,x)
}

BInv<-function(n,i,x){
  qbeta(x,i,n-i+1)
}

GInv<-function(n,i,x){
  qbeta(1-x,i,n-i+1)
}


B(20,2,0.1)
BInv(20,2,B(20,2,0.1))


G(20,2,0.1)
GInv(20,2,G(20,2,0.1))

n=10
t=0.6
l=sapply(1:n,function(i)BInv(n,i,t))
m=sapply(1:n,function(i)GInv(n,i,t))

l
m


plot(l)
lines(m)
