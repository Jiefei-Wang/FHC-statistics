%Compute the steck determinant using the algorithm from Miecznikowski(2017)
%Arguments:
%l: l_i's in the steck determinant
%u: m_i's in the steck determinant
%return:
%prob: steck's determinant
function prob=computeOrderProb(l,u)
N=size(l);
N=N(1);
firsRow=SMatrix_row_simple(l,u,1,1);
for i = 2:N
   firsRow(i:N)=firsRow(i:N)-firsRow(i-1)*SMatrix_row_simple(l,u,i,i);
   firsRow(i-1)=0;
end
prob=double(factorial(sym(N))*firsRow(N)*((mod(N,2)==1)*2-1));
