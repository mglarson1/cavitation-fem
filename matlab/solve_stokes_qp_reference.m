function [u,p,info] = solve_stokes_qp_reference(tri,x,y,dind,dval,tol,kkt_tol)
% Independent convex QP reference for full-gradient CR, mu=1, f=t=0.
% Minimize 1/2 u'*K*u subject to B*u>=0 and prescribed wall face means.
% quadprog multipliers for -B*u<=0 are the element pressures.
if nargin<6,tol=1e-12;end
if nargin<7,kkt_tol=1e-12;end
nt=size(tri,1);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[edges,~,ic]=unique(sort(E,2),'rows'); e=reshape(ic,nt,3); nu=2*size(edges,1);
% Direct affine-coordinate assembly, separately from the production solvers.
x1=x(tri(:,1));x2=x(tri(:,2));x3=x(tri(:,3));
y1=y(tri(:,1));y2=y(tri(:,2));y3=y(tri(:,3));
detJ=(x2-x1).*(y3-y1)-(x3-x1).*(y2-y1);area=abs(detJ)/2;
gx=-2*[y2-y3,y3-y1,y1-y2]./detJ;
gy=-2*[x3-x2,x1-x3,x2-x1]./detJ;
ix=zeros(nt,6);ix(:,1:2:6)=2*e-1;ix(:,2:2:6)=2*e;
ri=zeros(nt*18,1);ci=ri;vv=ri;k=0;
for a=1:3
 for b=1:3
  val=area.*(gx(:,a).*gx(:,b)+gy(:,a).*gy(:,b));
  for c=1:2
   ids=k+(1:nt);ri(ids)=ix(:,2*a-2+c);ci(ids)=ix(:,2*b-2+c);vv(ids)=val;k=k+nt;
  end
 end
end
K=sparse(ri,ci,vv,nu,nu);K=(K+K')/2;
bv=zeros(nt,6);bv(:,1:2:6)=area.*gx;bv(:,2:2:6)=area.*gy;
B=sparse(repmat((1:nt)',6,1),ix(:),bv(:),nt,nu);
[fixed,ia]=unique(dind(:));ud=dval(ia);free=setdiff((1:nu)',fixed);
H=K(free,free);g=K(free,fixed)*ud;Bf=B(:,free);bd=B(:,fixed)*ud;
opts=optimoptions('quadprog','Algorithm','interior-point-convex', ...
 'Display','off','OptimalityTolerance',tol,'ConstraintTolerance',tol, ...
 'StepTolerance',1e-15,'MaxIterations',200);
[uf,~,flag,out,lam]=quadprog(H,g,-Bf,bd,[],[],[],[],[],opts);
assert(flag>0,'QP did not converge, flag %d',flag)
u=zeros(nu,1);u(fixed)=ud;u(free)=uf;p=lam.ineqlin;div=(B*u)./area;
info=struct('exitflag',flag,'iterations',out.iterations,'K',K,'B',B, ...
 'area',area,'edges',edges,'free',free,'fixed',fixed,'div',div, ...
 'stationarity',norm(H*uf+g-Bf'*p,inf), ...
 'min_pressure',min(p),'min_divergence',min(div), ...
 'complementarity',max(abs(p.*div)),'energy',0.5*u'*K*u);
% Polish the QP estimate with equality-constrained working-set solves.
% Retain the current classification for near-zero constraints; switch only
% on definite sign violations. This prevents roundoff-driven cycling.
info.qp_u=u;info.qp_p=p;info.qp_div=div;
act=(100*p>=div);converged=false;
for it=1:100
 Ba=Bf(act,:);rhs=[-g;bd(act)];
 z=[H,-Ba';-Ba,sparse(nnz(act),nnz(act))]\rhs;
 u(free)=z(1:numel(free));p=zeros(nt,1);p(act)=z(numel(free)+1:end);
 div=(B*u)./area;
 next=act;next(p < -kkt_tol/100)=false;next(div < -kkt_tol)=true;
 if isequal(next,act),converged=true;break;end
 act=next;
end
info.polished_u=u;info.polished_p=p;info.polished_div=div;info.polished_act=act;
info.polished_iterations=it;info.polished_converged=converged;info.kkt_tolerance=kkt_tol;
info.polished_stationarity=norm(K(free,:)*u-Bf'*p,inf);
info.polished_min_pressure=min(p);info.polished_min_divergence=min(div);
info.polished_complementarity=max(abs(p.*div));
info.polished_energy=0.5*u'*K*u;
assert(converged && min(p)>=-kkt_tol/100 && min(div)>=-kkt_tol, ...
 'QP polishing failed the KKT sign checks.');
assert(info.polished_stationarity<1e-11 && info.polished_complementarity<1e-11, ...
 'QP polishing failed the stationarity/complementarity checks.');
end
