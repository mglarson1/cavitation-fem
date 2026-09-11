function table11 = reproduce_solver_table(outdir)
% Reproduce manuscript Table 11: the steepest pit with the full-gradient
% Crouzeix--Raviart pair, no jump penalty and traction-free ends.  The
% active-set iteration is stopped by the sign rule (tau = 1e-12); the
% smoothed Newton method uses s = 1e-8, undamped, and marks a sample
% cavitated when p_T < (s/gamma)^(1/2).

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

meshes=[128 8;256 16;512 32];n=size(meshes,1);gam=100;s=1e-8;
opts=struct('maxit',60,'tol',1e-10,'verbose',0,'damp',0);
[elements,active_iterations,active_peak,active_cavity,active_min_pressure, ...
    active_min_divergence,newton_iterations,newton_peak,newton_cavity]=deal(zeros(n,1));
for k=1:n
    [tri,x,z,dind,dval,Hf]=pit_mesh(meshes(k,1),meshes(k,2));
    elements(k)=size(tri,1);
    [~,~,~,p,iA]=evalc('solvedisccr(tri,x,z,gam,dind,dval,1e-12)');
    assert(iA.converged)
    m=metrics(tri,x,z,p,~iA.act,Hf);
    active_iterations(k)=iA.iterations;active_peak(k)=m.peak;active_cavity(k)=m.cavity;
    active_min_pressure(k)=min(p);active_min_divergence(k)=min(iA.div);
    [~,~,ps,iN]=solvedisccrn(tri,x,z,gam,s,dind,dval,opts);
    assert(iN.converged)
    m=metrics(tri,x,z,ps,ps<iN.cavthreshold,Hf);
    newton_iterations(k)=iN.iterations;newton_peak(k)=m.peak;newton_cavity(k)=m.cavity;
end
table11=table(elements,active_iterations,active_peak,active_cavity, ...
    active_min_pressure,active_min_divergence,newton_iterations,newton_peak,newton_cavity);
writetable(table11,fullfile(outdir,'solver_comparison.csv'));

check_close('Table 11 active-set iterations',active_iterations,[8;9;9],0,0);
check_close('Table 11 active-set peaks',active_peak,[1.98695;2.09338;2.12389],6e-6,0);
check_close('Table 11 active-set cavities',active_cavity,[5.937;5.937;5.889],2e-4,0);
check_close('Table 11 Newton iterations',newton_iterations,[11;12;12],0,0);
check_close('Table 11 Newton peaks',newton_peak,[1.98695;2.09338;2.12389],6e-6,0);
check_close('Table 11 Newton cavities',newton_cavity,[2.439;2.346;2.343],2e-4,0);
assert(min(active_min_pressure)>=-1e-12 && min(active_min_divergence)>=-1e-12)
fprintf('Reproduced Table 11 in %s\n',outdir)
end

function row=metrics(tri,x,z,p,cav,Hf)
L=24;c=1;xs=linspace(0,L,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
pp=nan(size(xs));pp(ok)=p(eid(ok));cc=false(size(xs));cc(ok)=cav(eid(ok));
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peak',max(pp(w)),'cavity',dx*nnz(cc&w));
end

function [tri,x,z,dind,dval,Hf]=pit_mesh(nx,nz)
L=24;xp=12;c=1;V=1;
Hf=@(x) c*(1+exp(-(x-xp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);x=xh;z=c-zh.*Hf(xh);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
sliding=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
shaped=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
