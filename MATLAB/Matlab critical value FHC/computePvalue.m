%Compute the p-value for an observed FHC statistic
%Arguments:
%N: number of p-values
%fhc: the FHC statistic
%return:
%p: The p-value for the observed fhc statistic
function p=computePvalue(N,fhc)
l=computeL(N,fhc);
u=ones(N,1);
p=1-computeOrderProb(l,u);
end

