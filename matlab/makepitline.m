% pitline.pdf: pressure along the centerline y = 0.5 for the two
% discretizations of the Reynolds model -- the conforming mixed one of
% Section 2.1 and the multiplier-free stabilized method -- with the pit
% profile shown scaled for reference.  Refinement level 5, 6144 elements.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

delta=1;r=0.35;xp=1.5;yp=0.5;nx=96;ny=32;gam0=1e-3;
rho  =@(x,y) ((x-xp).^2+(y-yp).^2)/r^2;
dfun =@(x,y) 1+delta*exp(-rho(x,y));
ffun =@(x,y) delta*exp(-rho(x,y)).*(2*(x-xp)/r^2);      % f = -dd/dx
dxfun=@(x,y) -delta*exp(-rho(x,y)).*(2*(x-xp)/r^2);
dyfun=@(x,y) -delta*exp(-rho(x,y)).*(2*(y-yp)/r^2);

[tri,xnod,ynod]=meshrect(3,1,nx,ny);
tol=1e-12;
onb=@(x,y) abs(x)<tol|abs(x-3)<tol|abs(y)<tol|abs(y-1)<tol;

% ---- conforming mixed, piecewise linear ----
dn=find(onb(xnod,ynod));
[P1,lam,i1]=solvereynolds(tri,xnod,ynod,dfun,ffun,100,dn);

% ---- multiplier-free stabilized, piecewise quadratic ----
TR=triangulation(tri,xnod,ynod);ed=TR.edges;
bn=double(onb(xnod,ynod));
[t6,x6,y6,b6]=sqtriKb(tri,xnod,ynod,ed,bn);
dn6=find(onb(x6,y6));
[P2,i2]=solvereynolds2(t6,x6,y6,dfun,dxfun,dyfun,gam0,dn6);
fprintf('mixed: %d its, max P %.6f | stabilized: %d its, max P %.6f, min P %.2e\n', ...
    i1.iterations,max(P1),i2.iterations,max(P2),min(P2));

% ---- centerline ----
c1=find(abs(ynod-0.5)<1e-9);[xa,ia]=sort(xnod(c1));pa=P1(c1(ia));
c2=find(abs(y6  -0.5)<1e-9);[xb,ib]=sort(x6(c2));  pb=P2(c2(ib));
dprof=dfun(xa,0.5);
sc=max(pa)/max(dprof-1);                    % pit profile scaled for reference

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 13 7],'Color','w')
ax=axes('Position',[0.115 0.145 0.855 0.83]);hold on,box on
hp=plot(xa,(dprof-1)*sc,'-','Color',[0.68 0.68 0.68],'LineWidth',3.0);
h1=plot(xa,pa,'k-','LineWidth',1.3);
h2=plot(xb,pb,'r--','LineWidth',1.1);
plot([0 3],[0 0],'k:','LineWidth',0.5)
xlim([0 3]),ylim([-0.004 0.048])
xlabel('{\itx}'),ylabel('{\itP_h}')
legend([h1 h2 hp],{'mixed, {\itk} = 1','stabilised, {\itk} = 2', ...
    'pit profile (scaled)'},'Location','NorthWest','Box','off','FontSize',8)
set(ax,'FontSize',8,'Layer','top')
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pitline.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pitline.pdf'));
