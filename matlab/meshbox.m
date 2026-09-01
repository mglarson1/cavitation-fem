function [tet,xnod,ynod,znod]=meshbox(Lx,Ly,Lz,nx,ny,nz)
% Structured tetrahedral mesh of the box [0,Lx]x[0,Ly]x[0,Lz].
%
% Each of the nx*ny*nz hexahedra is split into six tetrahedra by the Kuhn
% (Freudenthal) subdivision, all six sharing the main diagonal of the cell.
% That subdivision is conforming: every interior face carries the same
% diagonal seen from either side, so no hanging nodes appear.
%
% This is the three dimensional counterpart of the quadrilateral mesh plus
% split into triangles used in maintri.m.
%
%   out : tet          6*nx*ny*nz x 4 connectivity, positively oriented
%         xnod,..      node coordinates

[X,Y,Z]=ndgrid(linspace(0,Lx,nx+1),linspace(0,Ly,ny+1),linspace(0,Lz,nz+1));
xnod=X(:);ynod=Y(:);znod=Z(:);

[I,J,K]=ndgrid(1:nx,1:ny,1:nz);
I=I(:);J=J(:);K=K(:);
ncell=numel(I);

% corner b of a cell, b = bx + 2*by + 4*bz
C=zeros(ncell,8);
for b=0:7
    bx=bitand(b,1);by=bitand(bitshift(b,-1),1);bz=bitand(bitshift(b,-2),1);
    C(:,b+1)=(I+bx)+(J+by-1)*(nx+1)+(K+bz-1)*(nx+1)*(ny+1);
end

% the six Kuhn tetrahedra, in corner bit numbering
KUHN=[0,1,3,7;0,1,5,7;0,2,3,7;0,2,6,7;0,4,5,7;0,4,6,7];

tet=zeros(6*ncell,4);
for t=1:6
    tet((t-1)*ncell+(1:ncell),:)=C(:,KUHN(t,:)+1);
end

% make every tetrahedron positively oriented
a=[xnod(tet(:,2))-xnod(tet(:,1)),ynod(tet(:,2))-ynod(tet(:,1)),znod(tet(:,2))-znod(tet(:,1))];
b=[xnod(tet(:,3))-xnod(tet(:,1)),ynod(tet(:,3))-ynod(tet(:,1)),znod(tet(:,3))-znod(tet(:,1))];
c=[xnod(tet(:,4))-xnod(tet(:,1)),ynod(tet(:,4))-ynod(tet(:,1)),znod(tet(:,4))-znod(tet(:,1))];
dj=dot(cross(a,b,2),c,2);
flip=dj<0;
tet(flip,[3,4])=tet(flip,[4,3]);
