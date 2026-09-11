% stokes3dfine.pdf (Figure 5): Stokes flow with cavitation in three
% dimensions.  Left, the pressure on the boundary of the 36x12x12 mesh with
% the cavitated elements at the outflow face marked; right, the pressure
% along the line y = z = 0.5 on the 18x6x6 and 36x12x12 meshes against the
% Poiseuille value 6-2x.
%
% Each boundary face is drawn exactly once, either in the pressure colour or,
% if its element is cavitated, in red.  Drawing the cavitated faces a second
% time on top of the coloured ones fails in vector output, because the
% painters renderer cannot order coplanar faces and interleaves the two
% layers.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

gam=100;op=struct('gamma1',1,'lamfac',-2/3,'verbose',0);
xs=linspace(0.002,2.998,3000)';
% the line y = z = 1/2 runs along mesh edges, so it is sampled just off it
pts=[xs,(0.5+1e-3)*ones(size(xs)),(0.5+1e-3)*ones(size(xs))];

% ---- 18x6x6, direct solver: centreline only ----
[tet,x,y,z,dind,dval]=box3dcase(18,6,6);
[~,~,~,p,info]=solvedisccrs3(tet,x,y,z,gam,dind,dval,op);
tid=pointLocation(triangulation(tet,[x y z]),pts);ok=~isnan(tid);
pl18=nan(size(xs));pl18(ok)=p(tid(ok));

% ---- 36x12x12, iterated penalty solver ----
[tet,x,y,z,dind,dval]=box3dcase(36,12,12);nele=size(tet,1);
[~,~,~,p,info]=solvedisccrs3_ip(tet,x,y,z,gam,dind,dval,op);
fprintf('%d elements, %d dof, %d iterations, %d cavitated\n', ...
    nele,info.ndof,info.iterations,info.ncav);
tid=pointLocation(triangulation(tet,[x y z]),pts);ok=~isnan(tid);
pl36=nan(size(xs));pl36(ok)=p(tid(ok));

fu=info.faces;onb=info.onbnd;nface=size(fu,1);
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[~,~,ic]=unique(sort(F,2),'rows');clear F
kk=(1:4*nele)';sel=onb(ic);
elemOf=zeros(nface,1);elemOf(ic(sel))=mod(kk(sel)-1,nele)+1;
bf=find(onb);cav=~info.act(elemOf(bf));
bfa=bf(~cav);bfc=bf(cav);

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 18 6.5],'Color','w')

% ---------------- left: pressure on the boundary ----------------
ax1=axes('Position',[0.075 0.13 0.38 0.80]);hold on
patch('Faces',fu(bfa,:),'Vertices',[x y z],'FaceVertexCData',p(elemOf(bfa)), ...
      'FaceColor','flat','EdgeColor','none')
patch('Faces',fu(bfc,:),'Vertices',[x y z], ...
      'FaceColor',[0.85 0.15 0.15],'EdgeColor',[0.45 0 0],'LineWidth',0.1)
colormap(ax1,parula);clim(ax1,[0 max(p(elemOf(bf)))])
cb=colorbar('eastoutside');cb.Label.String='{\itp}';cb.FontSize=8;
cb.Position=[0.455 0.18 0.016 0.70];
view(128,20),axis equal,axis tight,box on,grid off
xlabel('{\itx}'),ylabel('{\ity}'),zlabel('{\itz}')
set(ax1,'FontSize',8,'BoxStyle','full')

% ---------------- right: the centreline against Poiseuille ----------------
ax2=axes('Position',[0.585 0.17 0.385 0.76]);hold on,box on
plot(xs,6-2*xs,'k--','LineWidth',1.0)
stairs(xs,pl18,'-','Color',[0.6 0.6 0.6],'LineWidth',0.9)
stairs(xs,pl36,'r-','LineWidth',1.1)
xlabel('{\itx}'),ylabel('{\itp}')
legend({'6-2{\itx}','18\times6\times6','36\times12\times12'}, ...
    'Location','NorthEast','Box','off','FontSize',8)
xlim([0 3]),ylim([-0.3 6.5]),set(ax2,'FontSize',8)

set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'stokes3dfine.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'stokes3dfine.pdf'));
