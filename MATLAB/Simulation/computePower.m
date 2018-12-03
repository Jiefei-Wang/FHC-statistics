function power=computePower(critical,N,F_alt,l_func)
u=ones(N,1);
l=F_alt(l_func(N,critical));
power=1-computeOrderProb(l,u);