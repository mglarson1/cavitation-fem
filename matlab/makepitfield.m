% pitfieldends.pdf (Figure 7): the Stokes solution near the steepest pit,
% delta/r = 1, with pressure-normal-flow ends (requires Optimization Toolbox).
% Above, the horizontal velocity with the contour u_x = 0 bounding the
% reversed flow in the pit; below, the pressure.  The cavitated elements are
% lightened in both.  The pit is a depression in the stationary surface
% below; the flat plane above slides to the right.
%
% Rendering note.  The velocity is Crouzeix-Raviart, so it is piecewise
% linear and discontinuous, with its degrees of freedom at the edge
% midpoints.  For the picture it is evaluated at the element vertices, where
% u|_T = sum_i u_i - 2 u_k at vertex k, and averaged over the elements
% meeting there; that gives a continuous field to contour.  The pressure is
% elementwise constant and is drawn flat, which is what it is.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

L=24;xp=L/2;c=1;V=1;gam=100;delta=1;rp=1;
nx=256;nz=16;
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));

[tri,xh,zh]=meshrect(L,1,nx,nz);
xm=xh; zm=c-zh.*Hf(xh);            % pit downwards, sliding plane at z = c

E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
flat=find(onb & abs(zh(eu(:,1)))<1e-9   & abs(zh(eu(:,2)))<1e-9);   % sliding
shap=find(onb & abs(zh(eu(:,1))-1)<1e-9 & abs(zh(eu(:,2))-1)<1e-9); % stationary
% pressure-normal-flow ends: transverse velocity prescribed, normal free
ends=find(onb & ((abs(xh(eu(:,1)))<1e-9 & abs(xh(eu(:,2)))<1e-9) | ...
    (abs(xh(eu(:,1))-L)<1e-9 & abs(xh(eu(:,2))-L)<1e-9)));
dind=[2*flat-1;2*flat;2*shap-1;2*shap;2*ends];
dval=[V*ones(numel(flat),1);zeros(numel(flat),1); ...
      zeros(numel(shap),1);zeros(numel(shap),1);zeros(numel(ends),1)];
dc=false(size(eu,1),2);dc([flat;shap],:)=true;dc(ends,2)=true;
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','qp', ...
    'dirichlet_components',dc,'kkt_tolerance',1e-11);
[ux,uy,p,info]=solvedisccrs(tri,xm,zm,gam,dind,dval,opCR);
info.ncav=nnz(~info.act);
fprintf('%d elements, p in [%.4g %.4g], %d cavitated\n', ...
    size(tri,1),min(p),max(p),info.ncav);

% ---- CR velocity at the vertices, averaged over the adjoining elements ----
edg=reshape(ic,size(tri,1),3);
nno=numel(xm);
acc=zeros(nno,1);cnt=zeros(nno,1);
for iel=1:size(tri,1)
    ue=ux(edg(iel,:));s=sum(ue);
    for k=1:3
        v=tri(iel,k);
        acc(v)=acc(v)+(s-2*ue(k));cnt(v)=cnt(v)+1;
    end
end
uxv=acc./max(cnt,1);

% the mesh is structured, x varying fastest, so reshape and contour directly
X=reshape(xm,nx+1,nz+1);Z=reshape(zm,nx+1,nz+1);UX=reshape(uxv,nx+1,nz+1);

x0=8;x1=16;
% only the elements inside the plotted window are drawn.  The rest would be
% clipped by the axes anyway, but in vector output they are still written to
% the file, and they are two thirds of the mesh.
xt=reshape(xm(tri),[],3);
vis=find(max(xt,[],2)>x0-0.2 & min(xt,[],2)<x1+0.2);
cav=intersect(find(~info.act),vis);
% the contour is likewise built on the visible columns only
col=find(X(:,1)>x0-0.2 & X(:,1)<x1+0.2);
Xc=X(col,:);Zc=Z(col,:);UXc=UX(col,:);
zlo=min(zm)-0.08;zhi=c+0.08;

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 17 10],'Color','w')
xx=linspace(0,L,1601)';
xtop=linspace(x0,x1,2)';
% Everything is drawn in a strictly two dimensional axes, later objects on
% top of earlier ones.  Giving the patches a z coordinate makes the axes
% three dimensional, and the painters renderer then decomposes the Gouraud
% shaded triangles into bands with white seams between them.

% ---------------- above: horizontal velocity ----------------
ax1=subplot(2,1,1);hold on
% Flat, one color per element, not interpolated.  MATLAB approximates
% Gouraud shading in vector output by a band decomposition that leaves white
% seams across the field, whichever object carries it; flat shading has no
% such problem and is in any case the honest rendering of a discontinuous
% piecewise linear velocity.  The element values are the means of the three
% vertex values, and the contour below still uses the vertex field.
uxe=mean(reshape(uxv(tri),[],3),2);
patch('Faces',tri(vis,:),'Vertices',[xm zm], ...
      'FaceVertexCData',uxe(vis),'FaceColor','flat','EdgeColor','none')
colormap(ax1,parula)
cb=colorbar('eastoutside');cb.Label.String='{\itu_x}';cb.FontSize=8;
patch('Faces',tri(cav,:),'Vertices',[xm zm], ...
      'FaceColor','w','FaceAlpha',0.55,'EdgeColor','none')
hc=contour(Xc,Zc,UXc,[0 0],'w','LineWidth',1.2);
% the walls last and on top: u_x vanishes on the stationary surface, so the
% zero contour runs along it, and the wall line covers that coincident part,
% leaving visible only the curve that encloses the reversed flow
plot(xtop,c*ones(2,1),'k-','LineWidth',1.0)
plot(xx,c-Hf(xx),'k-','LineWidth',1.0)
axis equal, axis([x0 x1 zlo zhi])
set(ax1,'FontSize',9,'Layer','top'),box on
ylabel('{\itz}'),set(ax1,'XTickLabel',[])
text(0.012,0.84,'{\itu_x}','Units','normalized','FontSize',9,'Color','k')

% ---------------- below: pressure ----------------
ax2=subplot(2,1,2);hold on
patch('Faces',tri(vis,:),'Vertices',[xm zm], ...
      'FaceVertexCData',p(vis),'FaceColor','flat','EdgeColor','none')
colormap(ax2,parula)
cb2=colorbar('eastoutside');cb2.Label.String='{\itp}';cb2.FontSize=8;
patch('Faces',tri(cav,:),'Vertices',[xm zm], ...
      'FaceColor','w','FaceAlpha',0.55,'EdgeColor','none')
plot(xtop,c*ones(2,1),'k-','LineWidth',1.0)
plot(xx,c-Hf(xx),'k-','LineWidth',1.0)
axis equal, axis([x0 x1 zlo zhi])
set(ax2,'FontSize',9,'Layer','top'),box on
xlabel('{\itx}'),ylabel('{\itz}')
text(0.012,0.84,'{\itp}','Units','normalized','FontSize',9,'Color','k')

% all text in the figure in the body font of the paper, which elsarticle
% takes from txfonts, so that axis labels and captions do not clash
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pitfieldends.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pitfieldends.pdf'));
