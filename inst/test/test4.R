
n=10000
stat=0.01

l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
m=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
mat = computeSMatrix(l,m)

mat2=computeLSMatrix(l,m)
range(mat2)/log(2)
system.time(
BJPvalue(stat,n,precBits = 512,autoPrecision=FALSE)
)
system.time(
  BJPvalue(stat,n,precBits = 1024*16,autoPrecision=FALSE)
)

BJCriticalMC(0.05,n)
BJPvalue(0.001132529,n)
