function [tri,xnod,ynod,quad]=meshrect(L,W,nx,ny)
% Structured triangulation of the rectangle (0,L) x (0,W) into nx by ny
% quadrilaterals, each split into two triangles.  Node numbering is
% column by column in x, so that xnod varies fastest.
%
%  out : tri    2*nx*ny x 3 triangles
%        quad   nx*ny x 4 quadrilaterals, counter clockwise
[XX,YY]=meshgrid(linspace(0,L,nx+1),linspace(0,W,ny+1));
xnod=reshape(XX',[],1);ynod=reshape(YY',[],1);   % x varies fastest
n=@(i,j) (j-1)*(nx+1)+i;                 % i = 1..nx+1 in x, j = 1..ny+1 in y
quad=zeros(nx*ny,4);k=0;
for j=1:ny
    for i=1:nx
        k=k+1;
        quad(k,:)=[n(i,j),n(i+1,j),n(i+1,j+1),n(i,j+1)];
    end
end
tri=[quad(:,1),quad(:,3),quad(:,4);quad(:,1),quad(:,2),quad(:,3)];
