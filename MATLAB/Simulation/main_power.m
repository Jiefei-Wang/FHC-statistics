clc;clear all;
pack;
addpath('symbolic',"method3",'boundFunc');

%Set the model parameter
%N: sample size, number of tests
N=100;
alpha=0.05;
accurate=0.0001;

%The alternative distribution beta(a,b)
a=0.5;
b=1;

%Compute the critical value for HC and FHC
[HC_critical,p]=gridSearch(N,alpha,@HC_computeL,1,100,accurate);
[FHC_critical,p]=gridSearch(N,alpha,@FHC_computeL,0,1,accurate);

%The proportion of alternative distribution
pi_list=0:0.01:1;
%The result will be stored in the following variable
record=zeros(length(pi_list),2);
for i=1:length(pi_list)
    i
    pi=pi_list(i);
    F_alt=@(x)pi*x+(1-pi)*betacdf(x,a,b);
    record(i,1)=computePower(FHC_critical,N,F_alt,@FHC_computeL);
    record(i,2)=computePower(HC_critical,N,F_alt,@HC_computeL);
end



