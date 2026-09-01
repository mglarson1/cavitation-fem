function verification = verify_newton_safeguard(outdir)
% Verify the optional residual-decrease safeguard on all 48 Newton cases.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end
reference_file=fullfile(outdir,'table10_smoothed_newton_all_48_cases.csv');
if(~exist(reference_file,'file'))
    error('verify_newton_safeguard:reference', ...
        'Run reproduce_smoothed_newton_table before this verification')
end
reference=readtable(reference_file);

meshes=[64 4;128 8;256 16;512 32];
svals=[1e-2;1e-4;1e-6;1e-8;1e-10;1e-12];
cases=[1 2;1 1];
nrun=size(meshes,1)*numel(svals)*size(cases,1);
delta=zeros(nrun,1);radius=delta;steepness=delta;mesh_nx=delta;
mesh_nz=delta;elements=delta;smoothing=delta;iterations=delta;
converged=false(nrun,1);backtracks=delta;pmax=delta;cavity=delta;
cavity_threshold=delta;central_path_error=delta;relative_residual=delta;

opts=struct('maxit',60,'tol',1e-10,'verbose',0,'damp',1);irun=0;
for icase=1:size(cases,1)
    for imesh=1:size(meshes,1)
        [tri,x,z,dind,dval,Hf]=pit_problem(cases(icase,1),cases(icase,2), ...
            meshes(imesh,1),meshes(imesh,2));
        for is=1:numel(svals)
            irun=irun+1;s=svals(is);
            [~,~,p,info]=solvedisccrn(tri,x,z,100,s,dind,dval,opts);
            metric=smoothed_metrics(tri,x,z,p,Hf,info.cavthreshold);
            delta(irun)=cases(icase,1);radius(irun)=cases(icase,2);
            steepness(irun)=delta(irun)/radius(irun);
            mesh_nx(irun)=meshes(imesh,1);mesh_nz(irun)=meshes(imesh,2);
            elements(irun)=size(tri,1);smoothing(irun)=s;
            iterations(irun)=info.iterations;converged(irun)=info.converged;
            backtracks(irun)=info.nbacktrack;pmax(irun)=metric.peak;
            cavity(irun)=metric.cavity;cavity_threshold(irun)=info.cavthreshold;
            central_path_error(irun)=info.central_path_error;
            relative_residual(irun)=info.relres;
            assert(info.converged)
        end
    end
end

assert(height(reference)==nrun)
assert(all(reference.delta==delta & reference.radius==radius & ...
    reference.mesh_nx==mesh_nx & reference.mesh_nz==mesh_nz & ...
    reference.smoothing==smoothing))
undamped_pmax=reference.pmax;undamped_cavity=reference.cavity;
peak_difference=abs(pmax-undamped_pmax);
cavity_difference=abs(cavity-undamped_cavity);
verification=table(delta,radius,steepness,mesh_nx,mesh_nz,elements, ...
    smoothing,iterations,converged,backtracks,pmax,undamped_pmax, ...
    peak_difference,cavity,undamped_cavity,cavity_difference, ...
    cavity_threshold,central_path_error,relative_residual);
writetable(verification, ...
    fullfile(outdir,'table10_safeguarded_newton_all_48_cases.csv'));

used=backtracks>0;
assert(nnz(used)==4)
assert(all(steepness(used)==0.5 & smoothing(used)==1e-4))
assert(all(backtracks(used)>=7 & backtracks(used)<=9))
assert(max(peak_difference)<5e-9)
assert(max(cavity_difference)<1e-12)
assert(max(central_path_error)<1e-6)

total_cases=nrun;full_step_cases=nnz(~used);backtracked_cases=nnz(used);
min_backtracks=min(backtracks(used));max_backtracks=max(backtracks(used));
max_peak_difference=max(peak_difference);
max_cavity_difference=max(cavity_difference);
summary=table(total_cases,full_step_cases,backtracked_cases,min_backtracks, ...
    max_backtracks,max_peak_difference,max_cavity_difference);
writetable(summary,fullfile(outdir,'table10_safeguard_summary.csv'));
fprintf('Verified the Newton safeguard on all %d cases in %s\n',nrun,outdir)
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

function row=smoothed_metrics(tri,x,z,p,Hf,threshold)
L=24;c=1;xs=linspace(0,L,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
ps=nan(size(xs));ps(ok)=p(eid(ok));w=xs>4*c&xs<L-4*c;
dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(ps(w)<threshold));
end
