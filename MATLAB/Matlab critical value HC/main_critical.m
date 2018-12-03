clc;clear all;
pack;
addpath('symbolic');

%For an observed FHC statistic, compute its p-value
N=10;
alpha=0.05
accurate=0.0001;
[hc_critical,p]=gridSearch(N,alpha,accurate)


%Compute the critical value for a given sample size
N=190;
hc=3.82835;
p_value=computePvalue(N,hc)