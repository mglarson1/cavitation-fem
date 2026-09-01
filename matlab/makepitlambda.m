% pitlambda.pdf: the multiplier lambda_h for the pit of Fig. pit, from the
% conforming mixed discretization of Section 2.1.  It vanishes where the film
% is intact and is supported on the cavitated region, against whose free
% boundary it concentrates.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

delta=1;r=0.35;xp=1.5;yp=0.5;nx=96;ny=32;
rho =@(x,y) ((x-xp).^2+(y-yp).^2)/r^2;
dfun=@(x,y) 1+delta*exp(-rho(x,y));
ffun=@(x,y) delta*exp(-rho(x,y)).*(2*(x-xp)/r^2);

[tri,xnod,ynod]=meshrect(3,1,nx,ny);
tol=1e-12;
dn=find(abs(xnod)<tol|abs(xnod-3)<tol|abs(ynod)<tol|abs(ynod-1)<tol);
[P,lam,info]=solvereynolds(tri,xnod,ynod,dfun,ffun,100,dn);
fprintf('lambda in [%.4g %.4g], %d cavitated nodes\n',min(lam),max(lam),info.ncav);

% filled contours rather than an interpolated patch, for the reason given in
% makepit.m: Gouraud shading leaves seams in vector output
nlev=160;
X=reshape(xnod,nx+1,ny+1);Y=reshape(ynod,nx+1,ny+1);
Lg=reshape(lam,nx+1,ny+1);

th=linspace(0,2*pi,400);cx=xp+r*cos(th);cy=yp+r*sin(th);

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 15 4.6],'Color','w')
ax=axes('Position',[0.075 0.20 0.79 0.74]);hold on
contourf(X,Y,Lg,nlev,'LineStyle','none')
plot(cx,cy,'k--','LineWidth',0.9)
colormap(ax,parula);cb=colorbar('eastoutside');cb.FontSize=8;
% no label on the color bar: MATLAB renders \lambda from a symbol font that
% it does not embed in the PDF, and the caption names the quantity anyway
axis equal,axis([0 3 0 1]),box on
set(ax,'FontSize',8,'Layer','top')
xlabel('{\itx}'),ylabel('{\ity}')
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pitlambda.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pitlambda.pdf'));
