function [tet,x,y,z,dind,dval]=box3dcase(nx,ny,nz)
% Fig 7 / Table 9 geometry and boundary data on the nx x ny x nz Kuhn mesh
% of (0,3)x(0,1)x(0,1), as in reproduce_stokes_3d_table.m: no slip on
% y = 0,1, parabolic inflow (face means of y - y^2) on x = 0, u_z = 0 on the
% symmetry planes z = 0,1, and a free outflow boundary x = 3.
[tet,x,y,z]=meshbox(3,1,1,nx,ny,nz);
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');
onb=(accumarray(ic,1)==1);
Y=y(fu);X=x(fu);Z=z(fu);tol=1e-9;
isy=onb&(all(abs(Y)<tol,2)|all(abs(Y-1)<tol,2));
isin=onb&all(abs(X)<tol,2);
isz=onb&(all(abs(Z)<tol,2)|all(abs(Z-1)<tol,2));
my1=sum(Y,2)/3;
my2=(sum(Y.^2,2)+Y(:,1).*Y(:,2)+Y(:,1).*Y(:,3)+Y(:,2).*Y(:,3))/6;
gin=my1-my2;
fy=find(isy);fi=find(isin);fz=find(isz&~isy&~isin);
dind=[3*fy-2;3*fy-1;3*fy;3*fi-2;3*fi-1;3*fi;3*fz];
dval=[zeros(3*numel(fy),1);gin(fi);zeros(2*numel(fi),1);zeros(numel(fz),1)];
end
