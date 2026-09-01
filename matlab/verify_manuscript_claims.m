function results = verify_manuscript_claims(outdir,mode)
% Verify quantitative claims stated in the manuscript prose.
%
%   verify_manuscript_claims(outdir,'quick') checks the inexpensive two-
%   dimensional claims.  The 'full' mode also runs the finest gentle-pit
%   comparison and the three-dimensional gamma-invariance check.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(nargin<2 || isempty(mode)), mode='full'; end
if(~exist(outdir,'dir')), mkdir(outdir); end
mode=lower(string(mode));
if(~any(mode==["quick","full"]))
    error('verify_manuscript_claims:mode','Unknown mode "%s"',mode)
end

%% Reynolds nodal method: residuals, gamma invariance, and an independent QP solve
delta=1;r=0.35;xp=1.5;yp=0.5;
rho=@(x,y) ((x-xp).^2+(y-yp).^2)/r^2;
dfun=@(x,y) 1+delta*exp(-rho(x,y));
dxfun=@(x,y) -delta*(2*(x-xp)/r^2).*exp(-rho(x,y));
dyfun=@(x,y) -delta*(2*(y-yp)/r^2).*exp(-rho(x,y));
ffun=@(x,y) -dxfun(x,y);

ref=4;nx=3*2^ref;ny=2^ref;
[tri,x,y]=meshrect(3,1,nx,ny);
evalc('[P,lam,info]=solvereynolds(tri,x,y,dfun,ffun,100);');
free=find(~info.isdir);
kkt_vector=info.K*P-info.F-info.m.*lam;
kkt_residual=norm(kkt_vector(free),inf);
load_norm=norm(info.F,inf);
complementarity=max(abs(P.*lam));

gammas=[1e-4;1e4];ng=numel(gammas);
gamma_relative_difference=zeros(ng,1);gamma_iterations=zeros(ng,1);
gamma_active_match=false(ng,1);
for k=1:ng
    evalc('[Pg,~,ig]=solvereynolds(tri,x,y,dfun,ffun,gammas(k));');
    gamma_relative_difference(k)=norm(Pg-P)/norm(P);
    gamma_iterations(k)=ig.iterations;
    gamma_active_match(k)=isequal(ig.act,info.act);
end
assert(max(gamma_relative_difference)<1e-11)
assert(all(gamma_iterations==info.iterations) && all(gamma_active_match))

[Pfree,pgs_sweeps,pgs_kkt_residual]=projected_gauss_seidel( ...
    info.K(free,free),info.F(free),1e-13,200000);
Ppgs=zeros(size(P));Ppgs(free)=Pfree;
pgs_inf_difference=norm(Ppgs-P,inf);
pgs_active_match=isequal(Ppgs(free)<=1e-11,info.act(free));
energy=0.5*P'*info.K*P-info.F'*P;
energy_pgs=0.5*Ppgs'*info.K*Ppgs-info.F'*Ppgs;
pgs_energy_relative_difference=abs(energy_pgs-energy)/max(abs(energy),eps);
assert(pgs_kkt_residual<2e-13 && pgs_inf_difference<2e-11)
assert(pgs_active_match && pgs_energy_relative_difference<2e-11)

% The location claim uses refinement level 5 and the x-projection of all
% active nodes.  Stating the refinement and projection avoids implying
% that the two-dimensional active set is literally an interval.
ref_fine=5;nx_fine=3*2^ref_fine;ny_fine=2^ref_fine;
[tri_fine,x_fine,y_fine]=meshrect(3,1,nx_fine,ny_fine);
evalc('[P_fine,~,i_fine]=solvereynolds(tri_fine,x_fine,y_fine,dfun,ffun,100);');
cavity_x_min=min(x_fine(i_fine.act));cavity_x_max=max(x_fine(i_fine.act));
[~,imax]=max(P_fine);peak_x=x_fine(imax);
hx=3/nx_fine;
assert(abs(cavity_x_min-0.03)<hx && abs(cavity_x_max-1.28)<hx)
assert(abs(peak_x-1.95)<hx)

% The stated unstable gamma0=1 example is refinement 4 (1536 elements).
TR=triangulation(tri,x,y);
bnod=double(abs(x)<1e-12|abs(x-3)<1e-12|abs(y)<1e-12|abs(y-1)<1e-12);
[t6,x6,y6,~]=sqtriKb(tri,x,y,TR.edges,bnod);
evalc('[Pbad,ibad]=solvereynolds2(t6,x6,y6,dfun,dxfun,dyfun,1);');
unstable_converged=ibad.converged;
unstable_iterations=ibad.iterations;
unstable_pmin=min(Pbad);unstable_pmax=max(Pbad);
assert(~unstable_converged && unstable_iterations==50)
assert(abs(unstable_pmin+0.12)<1e-2)

reynolds_checks=table(kkt_residual,load_norm,complementarity, ...
    max(gamma_relative_difference),info.iterations,all(gamma_active_match), ...
    pgs_sweeps,pgs_kkt_residual,pgs_inf_difference,pgs_active_match, ...
    pgs_energy_relative_difference,cavity_x_min,cavity_x_max,peak_x, ...
    unstable_converged,unstable_iterations,unstable_pmin,unstable_pmax, ...
    'VariableNames',{'kkt_residual','load_norm','complementarity', ...
    'max_gamma_relative_difference','gamma_iterations','gamma_active_match', ...
    'pgs_sweeps','pgs_kkt_residual','pgs_inf_difference','pgs_active_match', ...
    'pgs_energy_relative_difference','cavity_x_min','cavity_x_max','peak_x', ...
    'unstable_gamma0_converged','unstable_gamma0_iterations', ...
    'unstable_gamma0_pmin','unstable_gamma0_pmax'});
writetable(reynolds_checks,fullfile(outdir,'prose_reynolds_checks.csv'));

%% Two-dimensional Stokes gamma invariance
ref=4;nx=3*2^ref;ny=2^ref;
[tri,x,y]=meshrect(3,1,nx,ny);[dind,dval]=channel_cr_boundary(tri,x,y);
op=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
evalc('[~,~,p_lo,i_lo]=solvedisccrs(tri,x,y,1e-4,dind,dval,op);');
evalc('[~,~,p_hi,i_hi]=solvedisccrs(tri,x,y,1e4,dind,dval,op);');
pressure_relative_difference=norm(p_hi-p_lo)/norm(p_lo);
active_match=isequal(i_hi.act,i_lo.act);
iterations_match=i_hi.iterations==i_lo.iterations;
cavitated_elements=i_lo.ncav;
assert(pressure_relative_difference<1e-10 && active_match && iterations_match)
assert(cavitated_elements==29 && i_lo.iterations==4)
stokes_gamma_checks=table(pressure_relative_difference,active_match, ...
    iterations_match,cavitated_elements,i_lo.iterations, ...
    'VariableNames',{'pressure_relative_difference','active_match', ...
    'iterations_match','cavitated_elements','iterations'});
writetable(stokes_gamma_checks,fullfile(outdir,'prose_stokes_gamma_checks.csv'));

%% Jump-penalty sweep on the steep pit
gamma1=[0;0.1;1;2;100];n=numel(gamma1);
pmax=zeros(n,1);cavity=pmax;iterations=pmax;converged=false(n,1);
[tri,x,z,dind,dval,Hf]=pit_problem(1,1,256,16);
for k=1:n
    op=struct('gamma1',gamma1(k),'lamfac',-2/3,'verbose',0,'solver','as');
    evalc('[~,~,p,ip]=solvedisccrs(tri,x,z,100,dind,dval,op);');
    m=cr_metrics(tri,x,z,p,ip,Hf,0);
    pmax(k)=m.peak;cavity(k)=m.cavity;iterations(k)=ip.iterations;
    converged(k)=ip.converged;
end
assert(all(converged))
assert(abs(pmax(1)-1.50)<1e-3 && cavity(1)==0)
assert(max(abs(pmax(2:4)-pmax(3))/pmax(3))<0.02)
assert(abs(pmax(5)-3.50)<1e-2)
stabilization_sweep=table(gamma1,iterations,converged,pmax,cavity);
writetable(stabilization_sweep,fullfile(outdir,'prose_stabilization_sweep.csv'));

%% Failure mode of the unstabilized active-set iteration
evalc('[~,~,p_unstable,i_unstable]=solvedisccr(tri,x,z,100,dind,dval);');
tail=i_unstable.cavhistory(41:end);
tail_mean=mean(tail);tail_min=min(tail);tail_max=max(tail);
tail_range=tail_max-tail_min;
assert(~i_unstable.converged && i_unstable.iterations==100)
assert(abs(tail_mean-3200)<150 && tail_range<150)
iteration=(1:numel(i_unstable.cavhistory))';
cavitated_elements_history=i_unstable.cavhistory;
writetable(table(iteration,cavitated_elements_history), ...
    fullfile(outdir,'prose_unstabilized_active_history.csv'));
unstabilized_summary=table(i_unstable.iterations,i_unstable.converged, ...
    tail_mean,tail_min,tail_max,tail_range, ...
    'VariableNames',{'iterations','converged','tail_mean','tail_min', ...
    'tail_max','tail_range'});
writetable(unstabilized_summary, ...
    fullfile(outdir,'prose_unstabilized_active_summary.csv'));

%% Smoothed Newton claims on the gentle and steep pits
[tri_g,x_g,z_g,dind_g,dval_g,Hg]=pit_problem(1,2,128,8);
evalc('[~,~,p_active,i_active]=solvedisccr(tri_g,x_g,z_g,100,dind_g,dval_g);');
m_active=cr_metrics(tri_g,x_g,z_g,p_active,i_active,Hg,0);
opts=struct('maxit',60,'tol',1e-10,'verbose',0,'damp',0);
evalc('[~,~,p_smooth,i_smooth]=solvedisccrn(tri_g,x_g,z_g,100,1e-12,dind_g,dval_g,opts);');
m_smooth=cr_metrics(tri_g,x_g,z_g,p_smooth,i_smooth,Hg,i_smooth.cavthreshold);
pressure_relative_difference=norm(p_smooth-p_active)/norm(p_active);
assert(i_active.converged && i_smooth.converged)
assert(abs(m_active.peak-3.52576)<1e-5 && abs(m_smooth.peak-m_active.peak)<1e-5)
assert(abs(m_active.cavity-4.812)<2e-4 && m_smooth.cavity==m_active.cavity)
assert(pressure_relative_difference<1e-7)
newton_gentle=table(i_active.iterations,i_smooth.iterations,m_active.peak, ...
    m_smooth.peak,m_active.cavity,m_smooth.cavity, ...
    pressure_relative_difference,i_smooth.central_path_error, ...
    'VariableNames',{'active_iterations','newton_iterations','active_pmax', ...
    'newton_pmax','active_cavity','newton_cavity', ...
    'pressure_relative_difference','central_path_error'});
writetable(newton_gentle,fullfile(outdir,'prose_newton_gentle_pit.csv'));

[tri_n,x_n,z_n,dind_n,dval_n]=pit_problem(1,1,256,16);
evalc('[~,~,~,i_order]=solvedisccrn(tri_n,x_n,z_n,100,1e-12,dind_n,dval_n,opts);');
residual=i_order.rhist(:);observed_order=nan(size(residual));
for k=3:numel(residual)
    observed_order(k)=log(residual(k)/residual(k-1))/ ...
        log(residual(k-1)/residual(k-2));
end
newton_iteration=(0:numel(residual)-1)';
writetable(table(newton_iteration,residual,observed_order), ...
    fullfile(outdir,'prose_newton_order.csv'));
superlinear_phase=observed_order(5:8);
phase_min=min(superlinear_phase);phase_max=max(superlinear_phase);
maximum_observed_order=max(observed_order,[],'omitnan');
final_observed_order=observed_order(end);
assert(i_order.converged && phase_min>1.2 && phase_max<1.7)
newton_order_summary=table(phase_min,phase_max,maximum_observed_order, ...
    final_observed_order,i_order.iterations,i_order.relres, ...
    'VariableNames',{'phase_min','phase_max','maximum_observed_order', ...
    'final_observed_order','iterations','relative_residual'});
writetable(newton_order_summary, ...
    fullfile(outdir,'prose_newton_order_summary.csv'));

gentle_discretization=table;
stokes3d_gamma_checks=table;
if(mode=="full")
    %% Fine-mesh gentle-pit comparison quoted after Table 8
    [tri_gf,x_gf,z_gf,dind_gf,dval_gf,Hgf]=pit_problem(1,2,512,32);
    op=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
    evalc('[~,~,p_cr,i_cr]=solvedisccrs(tri_gf,x_gf,z_gf,100,dind_gf,dval_gf,op);');
    m_cr=cr_metrics(tri_gf,x_gf,z_gf,p_cr,i_cr,Hgf,0);
    m_th=taylor_hood_pit_metrics(1,2,512,32);
    assert(i_cr.converged && m_th.info.converged)
    assert(abs(m_cr.peak-m_th.peak)<1e-4)
    assert(abs(m_cr.cavity-2.06)<1e-2 && abs(m_th.cavity-1.89)<1e-2)
    method=["proposed";"taylor-hood"];
    pmax=[m_cr.peak;m_th.peak];cavity=[m_cr.cavity;m_th.cavity];
    iterations=[i_cr.iterations;m_th.info.iterations];
    gentle_discretization=table(method,iterations,pmax,cavity);
    writetable(gentle_discretization, ...
        fullfile(outdir,'prose_gentle_pit_discretization.csv'));

    %% Three-dimensional gamma invariance on the intermediate channel mesh
    [tet,x3,y3,z3]=meshbox(3,1,1,18,6,6);
    [dind3,dval3]=channel3d_boundary(tet,x3,y3,z3);
    op3=struct('gamma1',1,'lamfac',-2/3,'verbose',0);
    evalc('[~,~,~,p3_lo,i3_lo]=solvedisccrs3(tet,x3,y3,z3,1,dind3,dval3,op3);');
    evalc('[~,~,~,p3_hi,i3_hi]=solvedisccrs3(tet,x3,y3,z3,1e6,dind3,dval3,op3);');
    pressure_relative_difference=norm(p3_hi-p3_lo)/norm(p3_lo);
    active_match=isequal(i3_hi.act,i3_lo.act);
    iterations_match=i3_hi.iterations==i3_lo.iterations;
    cavitated_elements=i3_lo.ncav;
    assert(pressure_relative_difference<1e-10 && active_match && iterations_match)
    assert(cavitated_elements==71)
    stokes3d_gamma_checks=table(pressure_relative_difference,active_match, ...
        iterations_match,cavitated_elements,i3_lo.iterations, ...
        'VariableNames',{'pressure_relative_difference','active_match', ...
        'iterations_match','cavitated_elements','iterations'});
    writetable(stokes3d_gamma_checks, ...
        fullfile(outdir,'prose_stokes3d_gamma_checks.csv'));
end

results=struct('reynolds',reynolds_checks,'stokes_gamma',stokes_gamma_checks, ...
    'stabilization',stabilization_sweep,'unstabilized',unstabilized_summary, ...
    'newton_gentle',newton_gentle,'newton_order', ...
    table(newton_iteration,residual,observed_order), ...
    'newton_order_summary',newton_order_summary, ...
    'gentle_discretization',gentle_discretization, ...
    'stokes3d_gamma',stokes3d_gamma_checks);
fprintf('Verified manuscript prose claims in %s (%s mode)\n',outdir,mode)
end

function [x,sweeps,residual]=projected_gauss_seidel(A,b,tol,maxit)
n=numel(b);x=zeros(n,1);g=-b;d=diag(A);residual=inf;
for sweeps=1:maxit
    for i=1:n
        xn=max(0,x(i)-g(i)/d(i));dx=xn-x(i);
        if(dx~=0), x(i)=xn;g=g+dx*A(:,i); end
    end
    free=x>1e-12;
    residual=max([0;abs(g(free));max(-g(~free),0)]);
    if(residual<tol), return, end
end
error('verify_manuscript_claims:pgs','Projected Gauss--Seidel did not converge')
end

function [dind,dval]=channel_cr_boundary(tri,x,y)
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
ya=y(eu(:,1));yb=y(eu(:,2));xa=x(eu(:,1));xb=x(eu(:,2));
gmean=(ya+yb)/2-(ya.^2+ya.*yb+yb.^2)/3;
wall=find(onb&((abs(ya)<1e-12&abs(yb)<1e-12)| ...
    (abs(ya-1)<1e-12&abs(yb-1)<1e-12)));
inflow=find(onb&abs(xa)<1e-12&abs(xb)<1e-12);
dind=[2*wall-1;2*wall;2*inflow-1;2*inflow];
dval=[zeros(2*numel(wall),1);gmean(inflow);zeros(numel(inflow),1)];
end

function [tri,x,z,dind,dval,Hf]=pit_problem(delta,rp,nx,nz)
L=24;xp=12;c=1;V=1;
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);x=xh;z=c-zh.*Hf(xh);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
sliding=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
shaped=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
end

function row=cr_metrics(tri,x,z,p,info,Hf,threshold)
L=24;c=1;xs=linspace(0,L,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
ps=nan(size(xs));ps(ok)=p(eid(ok));w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
if(threshold>0)
    cav=false(size(xs));cav(ok)=ps(ok)<threshold;
else
    cav=false(size(xs));cav(ok)=~info.act(eid(ok));
end
row=struct('peak',max(ps(w)),'cavity',dx*nnz(cav(w)));
end

function row=taylor_hood_pit_metrics(delta,rp,nx,nz)
L=24;xp=12;c=1;V=1;
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);TR=triangulation(tri,xh,zh);
[t6,x6,zh6,b6]=sqtriKb(tri,xh,zh,TR.edges,ones(size(xh)));
z6=c-zh6.*Hf(x6);nnop=numel(xh);
sliding=find(abs(zh6)<1e-9);shaped=find(abs(zh6-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
opts=struct('lamfac',-2/3);
[~,~,p,info]=solvedisc(t6,x6,z6,nnop,b6,100,dind,dval,opts);
xs=linspace(0,L,8001)';points=[xs,c-0.5*Hf(xs)];
TRm=triangulation(tri,xh,c-zh.*Hf(xh));eid=pointLocation(TRm,points);
ok=~isnan(eid);ps=nan(size(xs));
bc=cartesianToBarycentric(TRm,eid(ok),points(ok,:));
ps(ok)=sum(bc.*p(tri(eid(ok),:)),2);
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(ps(w)<=1e-8),'info',info);
end

function [dind,dval]=channel3d_boundary(tet,x,y,z)
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');onb=(accumarray(ic,1)==1);
Y=y(fu);X=x(fu);Z=z(fu);tol=1e-9;
isy=onb&(all(abs(Y)<tol,2)|all(abs(Y-1)<tol,2));
isin=onb&all(abs(X)<tol,2);
isz=onb&(all(abs(Z)<tol,2)|all(abs(Z-1)<tol,2));
my1=sum(Y,2)/3;
my2=(sum(Y.^2,2)+Y(:,1).*Y(:,2)+Y(:,1).*Y(:,3)+Y(:,2).*Y(:,3))/6;
gin=my1-my2;
fy=find(isy);fi=find(isin);fz=find(isz&~isy&~isin);
dind=[3*fy-2;3*fy-1;3*fy;3*fi-2;3*fi-1;3*fi;3*fz];
dval=[zeros(3*numel(fy),1);gin(fi);zeros(2*numel(fi),1);zeros(numel(fz),1)];
end
