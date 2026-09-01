% channel.pdf: the two dimensional Stokes problem of Section 4.4, computed
% with the stabilized Crouzeix-Raviart element.  Above, the pressure and the
% velocity field; below, the outflow region enlarged with the cavitated
% elements marked.  The counterpart of the plotting block of maintri.m,
% brought up to the discretization of the paper: the pressure is now
% elementwise constant, so it is drawn flat on each triangle, and the
% velocity lives at the edge midpoints rather than at the nodes.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

gam=100;nx=48;ny=16;
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
[tri,xnod,ynod]=meshrect(3,1,nx,ny);

% inflow as the exact edge mean of y(1-y)
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
ya=ynod(eu(:,1));yb=ynod(eu(:,2));xa=xnod(eu(:,1));xb=xnod(eu(:,2));
gmean=(ya+yb)/2-(ya.^2+ya.*yb+yb.^2)/3;
wall=find(onb&((abs(ya)<1e-12&abs(yb)<1e-12)|(abs(ya-1)<1e-12&abs(yb-1)<1e-12)));
inflow=find(onb&abs(xa)<1e-12&abs(xb)<1e-12);
dind=[2*wall-1;2*wall;2*inflow-1;2*inflow];
dval=[zeros(2*numel(wall),1);gmean(inflow);zeros(numel(inflow),1)];

[ux,uy,p,info]=solvedisccrs(tri,xnod,ynod,gam,dind,dval,opCR);
fprintf('%d elements, %d dof, %d iterations, p in [%.4g %.4g], %d cavitated\n', ...
    size(tri,1),info.ndof,info.iterations,min(p),max(p),info.ncav);

xe=info.xe;ye=info.ye;

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 17 11],'Color','w')

% ---------------- above: pressure and velocity ----------------
ax1=subplot(2,1,1);hold on
patch('Faces',tri,'Vertices',[xnod ynod],'FaceVertexCData',p, ...
      'FaceColor','flat','EdgeColor','none')
colormap(ax1,parula)
cb=colorbar('eastoutside');cb.Label.String='{\itp}';cb.FontSize=8;
% velocity at a subsample of the edge midpoints
sel=find(mod(1:numel(xe),7)'==1);
quiver(xe(sel),ye(sel),ux(sel),uy(sel),1.1,'k','LineWidth',0.4,'MaxHeadSize',0.7)
plot([0 3 3 0 0],[0 0 1 1 0],'k-','LineWidth',0.6)
axis equal, axis([0 3 0 1])
set(ax1,'FontSize',9,'Layer','top'),box on
xlabel('{\itx}'),ylabel('{\ity}')

% ---------------- below: the outflow region enlarged ----------------
ax2=subplot(2,1,2);hold on
x0=2.2;
patch('Faces',tri,'Vertices',[xnod ynod],'FaceVertexCData',p, ...
      'FaceColor','flat','EdgeColor',[0.75 0.75 0.75],'LineWidth',0.2)
colormap(ax2,parula)
cav=find(~info.act);
patch('Faces',tri(cav,:),'Vertices',[xnod ynod], ...
      'FaceColor',[0.85 0.15 0.15],'EdgeColor',[0.5 0 0],'LineWidth',0.2)
plot([x0 3 3 x0],[0 0 1 1],'k-','LineWidth',0.6)
axis equal, axis([x0 3 0 1])
set(ax2,'FontSize',9,'Layer','top'),box on
xlabel('{\itx}'),ylabel('{\ity}')
text(0.03,0.93,sprintf('%d cavitated elements of %d',info.ncav,size(tri,1)), ...
     'Units','normalized','FontSize',8,'BackgroundColor','w','Margin',1)

% all text in the figure in the body font of the paper, which elsarticle
% takes from txfonts, so that axis labels and captions do not clash
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'channel.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'channel.pdf'));
