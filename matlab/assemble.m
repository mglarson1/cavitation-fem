function [row,col,val,up] = assemble(elemat,eqs,row,col,val,up)
% assemble(elemat) sets up assembly, elemat = size of matrix storage
%
% assemble(neq) set up storage;
%
% assemble(elemat,eqs,row,col,val,up) store matrices;
% elemat = square element matrix to be stored
% eqs = equation numbers
% row = row positions of storage
% col = column positions of storage
% val = stored element matrix
% up = last use place in storage
%
if(nargin==1)
    row = zeros(elemat, 1);
    col = row;
    val = row;
    up = 0;
    return
end
eqs=eqs(:);len = length(eqs);
X = eqs(:, ones(1, len));
Y = X';
nn = len * len;
lo = up + 1;
up = up + nn;
row(lo:up) = X(:);
col(lo:up) = Y(:);
val(lo:up) = elemat(:);
