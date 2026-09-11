function c=pressure_end_pit(delta,rp,left,right)
% Reynolds strip and Stokes cross section for the pit
%   H(x) = 1 + delta exp(-(x-12)^2/rp^2)
% on [left,right], with pressure-normal-flow Stokes ends (u_z = 0, zero
% normal traction), longitudinal spacing 24/256, 16 Stokes layers and four
% Reynolds strip layers.  The Stokes problem is solved by quadratic
% programming with working-set refinement.  Returns window diagnostics on
% 4 < x < 20 (sampling 0.000375, upper Stokes trace) and full-length
% mid-plane profiles for plotting.
if(nargin<3), left=0; right=24; end
hx=24/256;xp=12;dx=.003/8;nz=16;ny=4;
xx=(left:hx:right)';nx=numel(xx)-1;
H=@(x)1+delta*exp(-((x-xp)/rp).^2);
mdH=@(x)2*delta*(x-xp)/rp^2.*exp(-((x-xp)/rp).^2);     % f = -dH/dx

[triR,xi,yR]=meshrect(nx,4*hx,nx,ny);xR=xx(round(xi)+1);
evalc('[P,~,iR]=solvereynolds(triR,xR,yR,@(x,y)H(x),@(x,y)mdH(x),1,find(xR==left|xR==right));');
assert(iR.converged&&iR.signok);
mid=find(abs(yR-2*hx)<1e-12);[xr,ix]=sort(xR(mid));mid=mid(ix);
pr=6*P(mid);acR=iR.act(mid);

[tri,xi,zh]=meshrect(nx,1,nx,nz);x=xx(round(xi)+1);z=1-zh.*H(x);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[ed,~,ic]=unique(sort(E,2),'rows');bd=accumarray(ic,1)==1;
slide=find(bd&zh(ed(:,1))==0&zh(ed(:,2))==0);
wall=find(bd&zh(ed(:,1))==1&zh(ed(:,2))==1);
ends=find(bd&((x(ed(:,1))==left&x(ed(:,2))==left)| ...
    (x(ed(:,1))==right&x(ed(:,2))==right)));
di=[2*slide-1;2*slide;2*wall-1;2*wall;2*ends];
dv=[ones(numel(slide),1);zeros(numel(slide)+2*numel(wall)+numel(ends),1)];
dc=false(size(ed,1),2);dc([slide;wall],:)=true;dc(ends,2)=true;
op=struct('gamma1',1,'lamfac',-2/3,'solver','qp','verbose',0, ...
    'dirichlet_components',dc,'kkt_tolerance',1e-11);
[ux,uz,pS,iS]=solvedisccrs(tri,x,z,100,di,dv,op);

xs=(4+dx:dx:20-dx)';pRs=interp1(xr,pr,xs);ib=discretize(xs,xr);
cR=acR(ib)&acR(ib+1);
TR=triangulation(tri,x,z);iu=pointLocation(TR,[xs,1-.5*H(xs)+1e-10]);
pSs=pS(iu);cS=~iS.act(iu);
xf=linspace(0,24,8001)';
ef=pointLocation(TR,[xf,1-.5*H(xf)+1e-10]);okf=~isnan(ef);
pSf=nan(size(xf));pSf(okf)=pS(ef(okf));
cSf=false(size(xf));cSf(okf)=~iS.act(ef(okf));
pRf=interp1(xr,pr,xf);ibf=discretize(xf,xr);vf=~isnan(ibf);
cRf=false(size(xf));cRf(vf)=acR(ibf(vf))&acR(ibf(vf)+1);

c=struct('delta',delta,'r',rp,'peak_reynolds',max(pRs),'peak_stokes',max(pSs), ...
    'profile_difference_pct',100*max(abs(pSs-pRs))/max(pRs), ...
    'front_reynolds',front(xs,cR,xp),'front_stokes',front(xs,cS,xp), ...
    'left_censored',cR(1)&&cS(1),'min_pressure',min(pS),'min_divergence',min(iS.div), ...
    'complementarity',iS.compl,'stationarity',iS.stationarity,'elements',size(tri,1), ...
    'xf',xf,'pRf',pRf,'pSf',pSf,'cRf',cRf,'cSf',cSf, ...
    'tri',tri,'x',x,'z',z,'ux',ux,'uz',uz,'p',pS,'info',iS);
end

function f=front(xs,c,xp)
f=NaN;if(any(c)), f=xs(find(c,1,'last'))-xp; end
end
