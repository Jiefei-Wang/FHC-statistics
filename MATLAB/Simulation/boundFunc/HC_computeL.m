%The lower bounds of the joint uniform distribution
%Please see the paper to see how to compute them from an observed HC
%statistic
function l=HC_computeL(N,hc)
l=c_func((1:N)/N,hc/N^0.5)';