%Computing one row of steck determinant
%Arguments:
%l_sym: l_i's in the steck determinant in symbolic format
%m_sym: m_i's in the steck determinant in symbolic format
%rowIndex: The row index
%colIndex: The column index, the return value is S(rowIndex,colIndex:)
%Return:
%s: a row of steck matrix, s=S(rowIndex,colIndex:)
function s=SMatrix_row_simple(l_sym,m_sym,rowIndex,colIndex)

N=size(l_sym);
N=N(1);


numIndex=colIndex:N;


k=numIndex-rowIndex+1;
mMatrix=m_sym(rowIndex);
lMatrix=l_sym(numIndex)';
temMatrix=mMatrix-lMatrix;
temMatrix=sym(temMatrix);

s=temMatrix.^k./factorial(sym(k));






