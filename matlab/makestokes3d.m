% stokes3d.pdf: Stokes flow with cavitation in three dimensions on the
% 18x6x6 mesh, computed with the stabilized Crouzeix-Raviart element on
% tetrahedra.  Left, the pressure on the boundary with the cavitated
% elements at the outflow face marked; right, the pressure along the line
% y = z = 0.5 against the Poiseuille value 6-2x.
%
% The pressure is elementwise constant, so each boundary face is drawn flat
% in the value of the tetrahedron behind it.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

gam=100;nx=18;ny=6;nz=6;
op=struct('gamma1',1,'lamfac',-2/3,'verbose',0);
[tet,xnod,ynod,znod]=meshbox(3,1,1,nx,ny,nz);
nele=size(tet,1);

% ---- faces, in the same order the solver builds them ----
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');
nface=size(fu,1);
cnt=accumarray(ic,1);onb=(cnt==1);
elemOf=zeros(nface,1);
for k=1:4*nele
    e=ic(k);
    if(onb(e)), elemOf(e)=mod(k-1,nele)+1; end
end

% ---- boundary conditions, as in Section 4.7 ----
Y=ynod(fu);X=xnod(fu);Z=znod(fu);tol=1e-9;
isy =onb&(all(abs(Y)<tol,2)|all(abs(Y-1)<tol,2));     % no slip walls
isin=onb&all(abs(X)<tol,2);                            % inflow
isz =onb&(all(abs(Z)<tol,2)|all(abs(Z-1)<tol,2));      % symmetry planes
my1=sum(Y,2)/3;
my2=(sum(Y.^2,2)+Y(:,1).*Y(:,2)+Y(:,1).*Y(:,3)+Y(:,2).*Y(:,3))/6;
gin=my1-my2;                                           % face mean of y(1-y)
fy=find(isy);fi=find(isin);fz=find(isz&~isy&~isin);
dind=[3*fy-2;3*fy-1;3*fy; 3*fi-2;3*fi-1;3*fi; 3*fz];
dval=[zeros(3*numel(fy),1); gin(fi);zeros(2*numel(fi),1); zeros(numel(fz),1)];

[ux,uy,uz,p,info]=solvedisccrs3(tet,xnod,ynod,znod,gam,dind,dval,op);
fprintf('%d elements, %d dof, %d iterations, p in [%.4g %.4g], %d cavitated\n', ...
    nele,info.ndof,info.iterations,min(p),max(p),info.ncav);

bf=find(onb);
pb=p(elemOf(bf));
cavb=bf(~info.act(elemOf(bf)));

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 18 6.5],'Color','w')

% ---------------- left: pressure on the boundary ----------------
% viewed from the outflow end, so that the cavitated layer against x = 3 is
% the face nearest the reader
ax1=axes('Position',[0.075 0.13 0.38 0.80]);hold on
patch('Faces',fu(bf,:),'Vertices',[xnod ynod znod],'FaceVertexCData',pb, ...
      'FaceColor','flat','EdgeColor','none')
patch('Faces',fu(cavb,:),'Vertices',[xnod ynod znod], ...
      'FaceColor',[0.85 0.15 0.15],'EdgeColor',[0.45 0 0],'LineWidth',0.15)
colormap(ax1,parula)
cb=colorbar('eastoutside');cb.Label.String='{\itp}';cb.FontSize=8;
cb.Position=[0.455 0.18 0.016 0.70];
view(128,20),axis equal,axis tight,box on,grid off
xlabel('{\itx}'),ylabel('{\ity}'),zlabel('{\itz}')
set(ax1,'FontSize',8,'BoxStyle','full')

% ---------------- right: the centerline against Poiseuille ----------------
ax2=axes('Position',[0.585 0.17 0.385 0.76]);hold on,box on
TRt=triangulation(tet,[xnod ynod znod]);
xs=linspace(0.002,2.998,1500)';
% the line y = z = 1/2 runs along mesh edges, so it is sampled just off it
pts=[xs,(0.5+1e-3)*ones(size(xs)),(0.5+1e-3)*ones(size(xs))];
tid=pointLocation(TRt,pts);ok=~isnan(tid);
pl=nan(size(xs));pl(ok)=p(tid(ok));
plot(xs,6-2*xs,'k--','LineWidth',1.0)
stairs(xs,pl,'r-','LineWidth',1.1)
xlabel('{\itx}'),ylabel('{\itp}')
legend({'6-2{\itx}','computed'},'Location','NorthEast','Box','off','FontSize',8)
xlim([0 3]),ylim([-0.3 6.5])
set(ax2,'FontSize',8)

% all text in the figure in the body font of the paper, which elsarticle
% takes from txfonts, so that axis labels and captions do not clash
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'stokes3d.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'stokes3d.pdf'));
