%Compute the p-value for an observed HC statistic
%Arguments:
%N: number of p-values
%hc: the HC statistic
%return:
%p: The p-value for the observed fhc statistic
function p=computePvalue(N,hc)
u=ones(N,1);
l=c_func((1:N)/N,hc/N^0.5)';
p=1-computeOrderProb(l,u);
end

