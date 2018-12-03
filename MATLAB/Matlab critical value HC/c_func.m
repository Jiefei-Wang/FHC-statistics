%The C function to compute the l_i's in the steck matrix, see Miecznikowski(2017)
function c=c_func(x,a)
c=(x+(a^2-a*(a^2+4*(1-x).*x).^0.5)/2)/(1+a^2);