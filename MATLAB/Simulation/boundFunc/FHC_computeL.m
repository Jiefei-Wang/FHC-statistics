%The lower bounds of the joint uniform distribution
%Please see the paper to see how to compute them from an observed FHC
%statistic
function l=FHC_computeL(N,fhc)
l=ones(N,1);
for i=1:N
    l(i)=betainv(fhc,i,N+1-i);
end