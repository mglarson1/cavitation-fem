function results = verify_newton_tolerance(outdir)
% Focused small-s tolerance audit; leaves the Table 10 reference outputs intact.
if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end
[tri,x,z,dind,dval,Hf]=pit_mesh(1,1,256,16);
svals=[1e-10,1e-12]; tolerances=[1e-10,1e-12,1e-14];
rows=zeros(numel(svals)*numel(tolerances),13); k=0;
for s=svals
    pref=[]; cref=NaN;
    for tol=tolerances
        k=k+1;
        opts=struct('maxit',25,'tol',tol,'verbose',0,'damp',0);
        [~,~,p,info]=solvedisccrn(tri,x,z,100,s,dind,dval,opts);
        metric=smoothed_metrics(tri,x,z,p,Hf,info.cavthreshold);
        if(isempty(pref)),pref=p;cref=metric.cavity;end
        rows(k,:)=[s,tol,info.iterations,info.converged,info.relres, ...
            info.central_path_error,info.central_path_error/(s/100), ...
            min(p),min(info.div),metric.peak,metric.cavity, ...
            metric.cavity-cref,norm(p-pref)/norm(pref)];
        fprintf('s=%g tol=%g iterations=%d converged=%d cavity=%.6f scaled_product_error=%.3e\n', ...
            s,tol,info.iterations,info.converged,metric.cavity,info.central_path_error/(s/100));
        results=array2table(rows(1:k,:),'VariableNames', ...
            {'smoothing','tolerance','iterations','converged','relative_residual', ...
            'central_path_error','scaled_product_error','min_pressure','min_divergence', ...
            'pmax','cavity','cavity_change','relative_pressure_change'});
        writetable(results,fullfile(outdir,'newton_tolerance_sensitivity_2026-09-08.csv'));
    end
end
end

function [tri,x,z,dind,dval,Hf]=pit_mesh(delta,rp,nx,nz)
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

function row=smoothed_metrics(tri,x,z,p,Hf,threshold)
L=24;c=1;xs=linspace(0,L,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
ps=nan(size(xs));ps(ok)=p(eid(ok));w=xs>4*c&xs<L-4*c;
dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(ps(w)<threshold));
end
