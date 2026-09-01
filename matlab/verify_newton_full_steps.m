function verification = verify_newton_full_steps(outdir)
% Verify undamped Newton on the four cases that trigger the optional guard.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end
nx=[64;128;256;512];nz=[4;8;16;32];n=numel(nx);
elements=zeros(n,1);iterations=elements;converged=false(n,1);
backtracks=elements;pmax=elements;cavity=elements;relative_residual=elements;
central_path_error=elements;
opts=struct('maxit',60,'tol',1e-10,'verbose',0,'damp',0);
for k=1:n
    [tri,x,z,dind,dval,Hf]=pit_mesh(nx(k),nz(k));
    [~,~,p,info]=solvedisccrn(tri,x,z,100,1e-4,dind,dval,opts);
    elements(k)=size(tri,1);iterations(k)=info.iterations;
    converged(k)=info.converged;backtracks(k)=info.nbacktrack;
    [pmax(k),cavity(k)]=smoothed_metrics( ...
        tri,x,z,p,Hf,info.cavthreshold);
    relative_residual(k)=info.relres;
    central_path_error(k)=info.central_path_error;
end
verification=table(nx,nz,elements,iterations,converged,backtracks,pmax, ...
    cavity,relative_residual,central_path_error);
writetable(verification,fullfile(outdir,'table10_undamped_newton_check.csv'));
assert(all(converged) && all(backtracks==0))
fprintf('Verified full Newton steps on all four guarded cases in %s\n',outdir)
end

function [tri,x,z,dind,dval,Hf]=pit_mesh(nx,nz)
L=24;xp=12;c=1;V=1;delta=1;rp=2;
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

function [peak,cavity]=smoothed_metrics(tri,x,z,p,Hf,threshold)
xs=linspace(0,24,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,1-0.5*Hf(xs)]);ok=~isnan(eid);
ps=nan(size(xs));ps(ok)=p(eid(ok));w=xs>4&xs<20;
peak=max(ps(w));cavity=(xs(2)-xs(1))*nnz(ps(w)<threshold);
end
