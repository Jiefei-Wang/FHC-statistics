%This function compute the FHC statistic from the p-values
%Arguments:
%p:The p-values
%return:
%fhc: the FHC statistic
function fhc=computeFHC(p)

sp=sort(p);
fhcList=zeros(length(p),1);

for i=1:length(p)
    fhcList(i)=betacdf(sp(i),i,length(p)-i+1);
end

fhc=min(fhcList);