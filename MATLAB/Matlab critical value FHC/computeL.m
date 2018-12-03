%This function compute l_i in the steck's determinant
%Arguments:
%N: number of p-values
%fhc: the FHC statistic
%return:
%l: a sequence of l_i
function l=computeL(N,fhc)
l=ones(N,1);
for i=1:N
    l(i)=betainv(fhc,i,N+1-i);
end