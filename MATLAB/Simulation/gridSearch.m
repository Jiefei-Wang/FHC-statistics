%Compute the critical value for the HC/FHC statistic, it uses bisection method
%to find the solution
%Arguments:
%N:Number of (tests') p-values
%alpha: significance level
%l_func: the function to compute the lower bound of the joint uniform distribution
%for computing the p-values of HC/FHC statistic 
%lower/upper: Define the searching space for the bisection method
%accurate: The acceptable error for the critical value
%Return:
%critical: the critical value of FHC statistic
%prob: The p-value associated with the critical value
function [critical,prob]=gridSearch(N,alpha,l_func,lower,upper,accurate)
%check condition
p1=computePvalue(N,lower,l_func);
p2=computePvalue(N,upper,l_func);
direction=sign(p2-p1);
if((p1-alpha)*(p2-alpha)<=0)
    while(upper-lower>accurate)
        error=upper-lower;
        temp=(lower+upper)/2;
        prob=computePvalue(N,temp,l_func);
        if(direction*(prob-alpha)<0)
            lower=temp;
        else
            if(direction*(prob-alpha)>0)
                upper=temp;
            else
                break;
            end
        end
    end
    if(direction~=sign(p2-p1))
        'Warning:result may not be accurate,lower or upper bound is not valid'
        critical=-(lower+upper)/2;
        prob=alpha;
    else
        critical=(lower+upper)/2;
        prob=alpha;
    end
else
    'invalid upper or lower bound'
    prob=0;
    critical=0;
end




