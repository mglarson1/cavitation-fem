function [fi,fix,fiy,fiz,vol]=basis3(x,y,z,xc,yc,zc)
% Scalar evaluation of the linear basis functions on a tetrahedron.
% The three dimensional counterpart of basis.m.
%
%   xc,yc,zc   the four vertex coordinates
%   fi         the four barycentric coordinates at (x,y,z)
%   fix,..     their (constant) derivatives
%   vol        the volume of the tetrahedron

J=[xc(2)-xc(1),xc(3)-xc(1),xc(4)-xc(1);
   yc(2)-yc(1),yc(3)-yc(1),yc(4)-yc(1);
   zc(2)-zc(1),zc(3)-zc(1),zc(4)-zc(1)];
detj=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) ...
    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) ...
    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1));
Ji=inv(J);

xi=Ji*[x-xc(1);y-yc(1);z-zc(1)];
fi=[1-sum(xi);xi];

% row i of G is grad(lambda_i)
G=[-sum(Ji,1);Ji];
fix=G(:,1);fiy=G(:,2);fiz=G(:,3);

vol=abs(detj)/6;
