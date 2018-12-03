clc;clear all;
pack;
addpath('symbolic');


%For an observed FHC statistic, compute its p-value
N=190;
fhc=0.00067;
p_value=computePvalue(N,fhc)

%Compute the critical value for a given sample size
N=100;
alpha=0.05
accurate=0.0001;
[fhc_critical,p]=gridSearch(N,alpha,accurate)


