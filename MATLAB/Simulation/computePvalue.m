%Compute the p-value for an observed HC statistic
%Arguments:
%N: number of p-values
%hc: the FHC statistic
%return:
%p: The p-value for the observed fhc statistic
function p=computePvalue(N,hc,l_func)
l=l_func(N,hc);
u=ones(N,1);
p=1-computeOrderProb(l,u);
end

