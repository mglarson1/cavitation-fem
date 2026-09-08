function row=run_pressure_end_case(cfg,outdir)
% Matched pressure/normal-flow ends, fixed-window fronts and local refinement.
if ~exist(outdir,'dir'),mkdir(outdir);end
hx=24/256;xp=12;dx=.003/8;
if ~isfield(cfg,'factor'),cfg.factor=1;end
if ~isfield(cfg,'global_refine'),cfg.global_refine=false;end
if ~isfield(cfg,'kkt_tolerance'),cfg.kkt_tolerance=1e-11;end
if ~isfield(cfg,'shift'),cfg.shift=false;end
m=cfg.factor;nz=16*m;ny=4*m;
base=(cfg.left:hx:cfg.right)';xx=base(1);
for j=1:numel(base)-1
 n=1;if cfg.global_refine || (base(j)>=6 && base(j+1)<=18),n=m;end
 xx=[xx;base(j)+(1:n)'*hx/n]; %#ok<AGROW>
end
if cfg.shift
 xx=[xx(xx<6);6;(6+hx/(2*m):hx/m:18-hx/(2*m))';18;xx(xx>18)];
end
nx=numel(xx)-1;H=@(x)1+exp(-(x-xp).^2);
[triR,xi,yR]=meshrect(nx,4*hx,nx,ny);xR=xx(round(xi)+1);
[P,lam,iR]=solvereynolds(triR,xR,yR,@(x,y)H(x), ...
 @(x,y)2*(x-xp).*exp(-(x-xp).^2),1,find(xR==cfg.left|xR==cfg.right));
assert(iR.converged&&iR.signok);
mid=find(abs(yR-2*hx)<1e-12);[xr,ix]=sort(xR(mid));mid=mid(ix);pr=6*P(mid);acR=iR.act(mid);
[tri,xi,zh]=meshrect(nx,1,nx,nz);x=xx(round(xi)+1);z=1-zh.*H(x);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];[ed,~,ic]=unique(sort(E,2),'rows');bd=accumarray(ic,1)==1;
slide=find(bd&zh(ed(:,1))==0&zh(ed(:,2))==0);wall=find(bd&zh(ed(:,1))==1&zh(ed(:,2))==1);
ends=find(bd&((x(ed(:,1))==cfg.left&x(ed(:,2))==cfg.left)| ...
 (x(ed(:,1))==cfg.right&x(ed(:,2))==cfg.right)));
di=[2*slide-1;2*slide;2*wall-1;2*wall;2*ends];dv=[ones(numel(slide),1);zeros(numel(slide)+2*numel(wall)+numel(ends),1)];
dc=false(size(ed,1),2);dc([slide;wall],:)=true;dc(ends,2)=true;
op=struct('gamma1',1,'lamfac',-2/3,'solver','qp','verbose',0,'dirichlet_components',dc, ...
 'kkt_tolerance',cfg.kkt_tolerance);
tic;[ux,uz,pS,iS]=solvedisccrs(tri,x,z,100,di,dv,op);elapsed=toc;
xs=(4+dx:dx:20-dx)';pRs=interp1(xr,pr,xs);ib=discretize(xs,xr);cR=acR(ib)&acR(ib+1);
TR=triangulation(tri,x,z);iu=pointLocation(TR,[xs,1-.5*H(xs)+1e-10]);
il=pointLocation(TR,[xs,1-.5*H(xs)-1e-10]);
pSs=pS(iu);cS=~iS.act(iu);cSl=~iS.act(il);
rr=front_metrics(xs,cR,dx,xp);ss=front_metrics(xs,cS,dx,xp);sl=front_metrics(xs,cSl,dx,xp);
row=table(cfg.left,cfg.right,m,double(cfg.global_refine),cfg.kkt_tolerance,nx,nz,ny, ...
 size(tri,1),max(pRs),max(pSs),rr(1),ss(1),rr(2),ss(2),rr(3),ss(3), ...
 rr(4),ss(4),rr(5),ss(5),ss(6),sl(2), ...
 iR.iterations,iS.qp.qp_iterations,iS.iterations,iS.stationarity,min(pS),min(iS.div),iS.compl,elapsed, ...
 'VariableNames',{'left','right','factor','global_refine','kkt_tolerance','nx','nz','ny_reynolds', ...
 'stokes_elements','reynolds_peak','stokes_peak','reynolds_length','stokes_length', ...
 'reynolds_front','stokes_front','reynolds_front_upper','stokes_front_upper', ...
 'reynolds_left_censored','stokes_left_censored','reynolds_segments','stokes_segments', ...
 'stokes_first_relative_x','stokes_lower_trace_front','reynolds_iterations', ...
 'qp_iterations','polish_iterations','stationarity','min_pressure','min_divergence','complementarity','stokes_seconds'});
profiles=table(xs,pRs,pSs,cR,cS,cSl,iS.div(iu),'VariableNames', ...
 {'x','reynolds_pressure','stokes_pressure','reynolds_cavity','stokes_cavity','stokes_lower_cavity','stokes_divergence'});
tag=sprintf('case_%g_%g_m%d_g%d_tol%g',cfg.left,cfg.right,m,cfg.global_refine,cfg.kkt_tolerance);
if cfg.shift,tag=[tag,'_shift'];end
row.shift=cfg.shift;
writetable(row,fullfile(outdir,[tag,'.csv']));writetable(profiles,fullfile(outdir,[tag,'_profiles.csv']));
save(fullfile(outdir,[tag,'.mat']),'cfg','row','profiles','tri','x','z','ux','uz','pS','iS', ...
 'triR','xR','yR','P','lam','iR','xx');
fprintf('MATCH [%g,%g] m%d global%d: elements%d peaks %.8f/%.8f lengths %.6f/%.6f fronts %.6f/%.6f minD%.2e %.1fs\n', ...
 cfg.left,cfg.right,m,cfg.global_refine,size(tri,1),max(pRs),max(pSs),rr(1),ss(1),rr(2),ss(2),min(iS.div),elapsed);
end

function a=front_metrics(xs,c,dx,xp)
ii=find(c);f=NaN;hi=NaN;first=NaN;
if ~isempty(ii)
 f=xs(ii(end))-xp;first=xs(ii(1))-xp;
 if ii(end)<numel(xs),hi=xs(ii(end)+1)-xp;end
end
a=[dx*nnz(c),f,hi,c(1),nnz(diff([false;c])==1),first];
end
