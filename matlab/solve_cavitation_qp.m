function [u,p,act,info]=solve_cavitation_qp(K,B,area,ind,bcval,force,opts)
% Convex velocity QP followed by tolerance-aware equality-constrained solves.
% KKT conditions: K*u-force-B'*p=0 (free dofs), p>=0, B*u/area>=0.
if ~isfield(opts,'qp_tolerance'),opts.qp_tolerance=1e-12;end
if ~isfield(opts,'kkt_tolerance'),opts.kkt_tolerance=1e-11;end
% Remove assembly roundoff asymmetry before forming the convex energy.
K=(K+K')/2;
nu=size(K,1);nt=size(B,1);free=setdiff((1:nu)',ind);H=K(free,free);
g=K(free,ind)*bcval-force(free);Bf=B(:,free);bd=B(:,ind)*bcval;
opt=optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
 'OptimalityTolerance',opts.qp_tolerance,'ConstraintTolerance',opts.qp_tolerance, ...
 'StepTolerance',1e-15,'MaxIterations',200);
[uf,~,flag,out,lm]=quadprog(H,g,-Bf,bd,[],[],[],[],[],opt);
assert(flag>0,'Interior-point QP did not converge');
u=zeros(nu,1);u(ind)=bcval;u(free)=uf;p=lm.ineqlin;D=B*u./area;
act=100*p>=D;converged=false;
for it=1:100
 Ba=Bf(act,:);z=[H,-Ba';-Ba,sparse(nnz(act),nnz(act))]\[-g;bd(act)];
 u(free)=z(1:numel(free));p=zeros(nt,1);p(act)=z(numel(free)+1:end);D=B*u./area;
 next=act;next(p < -opts.kkt_tolerance/100)=false;next(D < -opts.kkt_tolerance)=true;
 if isequal(next,act),converged=true;break;end
 act=next;
end
station=norm(K(free,:)*u-force(free)-Bf'*p,inf);
info=struct('qp_flag',flag,'qp_iterations',out.iterations,'iterations',it, ...
 'converged',converged,'stationarity',station,'min_pressure',min(p), ...
 'min_divergence',min(D),'complementarity',max(abs(p.*D)), ...
 'kkt_tolerance',opts.kkt_tolerance);
assert(converged && min(p)>=-opts.kkt_tolerance/100 && min(D)>=-opts.kkt_tolerance, ...
 'Working-set refinement did not pass the sign checks');
assert(station<1e-10 && max(abs(p.*D))<1e-10,'KKT residual check failed');
end
