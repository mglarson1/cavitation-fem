% pit.pdf: the pit of Section 4.1.  Film thickness d above, scaled
% pressure P_h below, with the broken circle of radius r about the center of
% the pit.  Computed with the conforming mixed discretization of Section 2.1
% on the refinement level 5 mesh, 6144 elements.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

delta=1;r=0.35;xp=1.5;yp=0.5;
nx=96;ny=32;
rho =@(x,y) ((x-xp).^2+(y-yp).^2)/r^2;
dfun=@(x,y) 1+delta*exp(-rho(x,y));
ffun=@(x,y) delta*exp(-rho(x,y)).*(2*(x-xp)/r^2);   % f = -dd/dx

[tri,xnod,ynod]=meshrect(3,1,nx,ny);
tol=1e-12;
dn=find(abs(xnod)<tol|abs(xnod-3)<tol|abs(ynod)<tol|abs(ynod-1)<tol);
[P,lam,info]=solvereynolds(tri,xnod,ynod,dfun,ffun,100,dn);
d=dfun(xnod,ynod);
% The mesh is structured, x varying fastest, so the nodal fields reshape to a
% grid and can be drawn as filled contours.  That is done instead of an
% interpolated patch because MATLAB approximates Gouraud shading in vector
% output by a band decomposition that leaves white seams across the field;
% filled contours are flat polygons throughout and vectorize cleanly, while
% at this many levels they are visually continuous.
nlev=160;
X=reshape(xnod,nx+1,ny+1);Y=reshape(ynod,nx+1,ny+1);
D=reshape(d,nx+1,ny+1);Pg=reshape(P,nx+1,ny+1);

th=linspace(0,2*pi,400);cx=xp+r*cos(th);cy=yp+r*sin(th);

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 15 8],'Color','w')

ax1=subplot(2,1,1);hold on
contourf(X,Y,D,nlev,'LineStyle','none')
plot(cx,cy,'k--','LineWidth',0.9)
colormap(ax1,parula);cb=colorbar('eastoutside');
cb.Label.String='{\itd}';cb.FontSize=8;
axis equal,axis([0 3 0 1]),box on
set(ax1,'FontSize',8,'Layer','top','XTickLabel',[])
ylabel('{\ity}')
text(0.012,0.86,'{\itd}','Units','normalized','FontSize',9)

ax2=subplot(2,1,2);hold on
contourf(X,Y,Pg,nlev,'LineStyle','none')
plot(cx,cy,'k--','LineWidth',0.9)
colormap(ax2,parula);cb2=colorbar('eastoutside');
cb2.Label.String='{\itP_h}';cb2.FontSize=8;
axis equal,axis([0 3 0 1]),box on
set(ax2,'FontSize',8,'Layer','top')
xlabel('{\itx}'),ylabel('{\ity}')
text(0.012,0.86,'{\itP_h}','Units','normalized','FontSize',9)

set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pit.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pit.pdf'));
