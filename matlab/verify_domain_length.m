function results=verify_domain_length(outdir,domains)
% End-location sensitivity at fixed pit geometry, mesh spacing and window.
% Steep pit from Table 6; proposed deviatoric, jump-stabilized CR Stokes.
% Bounds default to symmetric extensions of the original [0,24] domain.
if nargin<1,outdir=fullfile(pwd,'results');end
if nargin<2,domains=[0,24;-12,36;-36,60];end
if ~exist(outdir,'dir'),mkdir(outdir);end
hx=24/256;nz=16;ny=4;xp=12;delta=1;rp=1;dx=.003;
Hf=@(x) 1+delta*exp(-((x-xp)/rp).^2);
dfun=@(x,y) Hf(x);
ffun=@(x,y) 2*(x-xp)/rp^2.*delta.*exp(-((x-xp)/rp).^2);
xs=linspace(0,24,8001)';win=xs>4&xs<20;
opts=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
rows=[];reference=[];
for k=1:size(domains,1)
 left=domains(k,1);right=domains(k,2);L=right-left;nx=round(L/hx);
 assert(abs(nx*hx-L)<1e-12 && left<=0 && right>=24);
 tic;
 [triR,xR,yR]=meshrect(L,ny*hx,nx,ny);xR=xR+left;
 % The central longitudinal coordinates, strip width, triangle split and
 % mapped Stokes coordinates are identical to the original discretization.
 assert(isequal(unique(xR(xR>=0&xR<=24)),(0:256)'*hx));
 dnR=find(abs(xR-left)<1e-9|abs(xR-right)<1e-9);
 [P,lambdaR,iR]=solvereynolds(triR,xR,yR,dfun,ffun,1,dnR);
 pR=6*P;mid=find(abs(yR-ny*hx/2)<1e-9);
 [xr,ix]=sort(xR(mid));mid=mid(ix);pr=pR(mid);ar=iR.act(mid);
 [pRs,cR]=reynolds_samples(xr,pr,ar,xs);
 [tri,x,zh]=meshrect(L,1,nx,nz);x=x+left;z=1-zh.*Hf(x);
 E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
 [eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
 slide=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
 wall=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
 di=[2*slide-1;2*slide;2*wall-1;2*wall];
 dv=[ones(numel(slide),1);zeros(numel(slide)+2*numel(wall),1)];
 [ux,uz,pS,iS]=solvedisccrs(tri,x,z,100,di,dv,opts);
 TR=triangulation(tri,x,z);ei=pointLocation(TR,[xs,1-.5*Hf(xs)]);
 assert(all(~isnan(ei)));
 pSs=pS(ei);cS=~iS.act(ei);
 if isempty(reference),reference=[pRs,pSs,double(cR),double(cS)];end
 mr=cavity_summary(xs(win),cR(win),dx);ms=cavity_summary(xs(win),cS(win),dx);
 % Secondary diagnostic with moving end exclusions. Kept distinct from the
 % fixed window, since increasing the observation window changes its length.
 xa=(left:dx:right)';wa=xa>left+4&xa<right-4;
 [~,caR]=reynolds_samples(xr,pr,ar,xa);
 ia=pointLocation(TR,[xa,1-.5*Hf(xa)]);ok=~isnan(ia);caS=false(size(xa));caS(ok)=~iS.act(ia(ok));
 rfree=~iR.isdir;
 kr=norm(iR.K(rfree,:)*P-iR.F(rfree)-iR.m(rfree).*lambdaR(rfree),inf);
 rows=[rows;left,right,L,nx,nz,size(tri,1),iR.converged,iS.converged, ...
  iR.signok,iS.signok,iR.iterations,iS.iterations,max(pRs(win)),max(pSs(win)), ...
  mr,ms,dx*nnz(caR(wa)),dx*nnz(caS(wa)), ...
  100*norm(pRs(win)-reference(win,1),inf)/max(reference(win,1)), ...
  100*norm(pSs(win)-reference(win,2),inf)/max(reference(win,2)), ...
  dx*nnz(xor(cR(win),logical(reference(win,3)))), ...
  dx*nnz(xor(cS(win),logical(reference(win,4)))),min(P),min(lambdaR),kr, ...
  iR.complementarity,min(pS),min(iS.div),iS.compl,toc]; %#ok<AGROW>
 names={'left','right','length','nx','nz','stokes_elements','reynolds_converged', ...
  'stokes_converged','reynolds_signok','stokes_signok','reynolds_iterations', ...
  'stokes_iterations','reynolds_peak','stokes_peak', ...
  'reynolds_cavity','reynolds_first_cavity_x','reynolds_last_cavity_x', ...
  'reynolds_cavity_segments','reynolds_left_censored','reynolds_right_censored', ...
  'stokes_cavity','stokes_first_cavity_x','stokes_last_cavity_x', ...
  'stokes_cavity_segments','stokes_left_censored','stokes_right_censored', ...
  'reynolds_moving_window_cavity','stokes_moving_window_cavity', ...
  'reynolds_profile_change_pct','stokes_profile_change_pct', ...
  'reynolds_set_mismatch_length','stokes_set_mismatch_length', ...
  'reynolds_min_pressure','reynolds_min_multiplier','reynolds_equilibrium', ...
  'reynolds_complementarity','stokes_min_pressure','stokes_min_divergence', ...
  'stokes_complementarity','seconds'};
 results=array2table(rows,'VariableNames',names);
 writetable(results,fullfile(outdir,'domain_length_sensitivity.csv'));
 tag=sprintf('domain_%g_%g',left,right);
 profiles=table(xs,pRs,pSs,cR,cS,'VariableNames', ...
  {'x','reynolds_pressure','stokes_pressure','reynolds_cavitated','stokes_cavitated'});
 writetable(profiles,fullfile(outdir,[tag,'_profiles.csv']));
 save(fullfile(outdir,[tag,'.mat']),'tri','x','z','ux','uz','pS','iS', ...
  'triR','xR','yR','P','lambdaR','iR','profiles','left','right');
 fprintf('Domain [%g,%g] %d elements: convergence R%d S%d; peak R%.8f S%.8f; cavity R%.6f S%.6f; clipped R%d%d S%d%d; %.1fs\n', ...
  left,right,size(tri,1),iR.converged,iS.converged,max(pRs(win)),max(pSs(win)), ...
  mr(1),ms(1),mr(5),mr(6),ms(5),ms(6),rows(end,end));
end
end

function [p,c]=reynolds_samples(xr,pr,ar,xs)
p=interp1(xr,pr,xs);ib=discretize(xs,xr);ok=~isnan(ib);c=false(size(xs));
c(ok)=ar(ib(ok))&ar(ib(ok)+1);
end

function v=cavity_summary(xs,c,dx)
idx=find(c);first=NaN;last=NaN;
if ~isempty(idx),first=xs(idx(1));last=xs(idx(end));end
v=[dx*nnz(c),first,last,nnz(diff([false;c])==1),c(1),c(end)];
end
