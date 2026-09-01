function results = reproduce_smoothed_newton_table(outdir)
% Reproduce manuscript Table 10 and all 48 smoothed-Newton experiments.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

meshes=[64 4;128 8;256 16;512 32];
svals=[1e-2;1e-4;1e-6;1e-8;1e-10;1e-12];
cases=[1 2;1 1];                         % delta/r = 0.5 and 1
nrun=size(meshes,1)*numel(svals)*size(cases,1);
delta=zeros(nrun,1);radius=delta;steepness=delta;mesh_nx=delta;
mesh_nz=delta;elements=delta;smoothing=delta;iterations=delta;
converged=false(nrun,1);backtracks=delta;pmax=delta;cavity=delta;
cavity_threshold=delta;central_path_error=delta;relative_residual=delta;

% Use the undamped Newton method reported in the manuscript.  The solver's
% optional residual-decrease safeguard is tested separately by callers.
opts=struct('maxit',60,'tol',1e-10,'verbose',0,'damp',0);irun=0;
for icase=1:size(cases,1)
    for imesh=1:size(meshes,1)
        [tri,x,z,dind,dval,Hf]=pit_mesh(cases(icase,1),cases(icase,2), ...
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
all48=table(delta,radius,steepness,mesh_nx,mesh_nz,elements,smoothing, ...
    iterations,converged,backtracks,pmax,cavity,cavity_threshold, ...
    central_path_error,relative_residual);
writetable(all48,fullfile(outdir,'table10_smoothed_newton_all_48_cases.csv'));
assert(height(all48)==48 && all(all48.converged) && all(all48.backtracks==0))
assert(max(all48.central_path_error)<1e-6)

%% Left half of Table 10: mesh dependence at s=1e-8, steep pit
sel=all48.steepness==1&all48.smoothing==1e-8;
smooth=all48(sel,:);smooth=sortrows(smooth,'elements');
active_iterations=zeros(4,1);active_converged=false(4,1);
for k=1:4
    [tri,x,z,dind,dval]=pit_mesh(1,1,meshes(k,1),meshes(k,2));
    evalc('[~,~,~,iA]=solvedisccr(tri,x,z,100,dind,dval);');
    active_iterations(k)=iA.iterations;active_converged(k)=iA.converged;
end
table10_mesh=table(smooth.mesh_nx,smooth.mesh_nz,smooth.elements, ...
    active_iterations,active_converged,smooth.iterations,smooth.pmax, ...
    smooth.cavity,'VariableNames',{'mesh_nx','mesh_nz','elements', ...
    'active_set_iterations','active_set_converged','newton_iterations', ...
    'pmax','cavity'});
writetable(table10_mesh,fullfile(outdir,'table10_smoothed_newton_mesh.csv'));
check_close('Table 10 active iterations',active_iterations,[8;9;100;100],0,0);
assert(all(active_converged(1:2)) && all(~active_converged(3:4)))
check_close('Table 10 Newton iterations by mesh',smooth.iterations,[10;11;12;12],0,0);
check_close('Table 10 peaks by mesh',smooth.pmax, ...
    [1.7087;1.9870;2.0934;2.1239],6e-5,0);

%% Right half of Table 10: smoothing dependence on 8192 elements
sel=all48.steepness==1&all48.elements==8192&all48.smoothing<=1e-4;
table10_smoothing=sortrows(all48(sel,:), 'smoothing','descend');
writetable(table10_smoothing, ...
    fullfile(outdir,'table10_smoothed_newton_smoothing.csv'));
check_close('Table 10 Newton iterations by smoothing', ...
    table10_smoothing.iterations,[15;14;12;10;10],0,0);
check_close('Table 10 peaks by smoothing',table10_smoothing.pmax, ...
    [2.09376;2.09339;2.09338;2.09338;2.09338],6e-6,0);
check_close('Table 10 cavities by smoothing',table10_smoothing.cavity, ...
    [0.753;1.689;2.346;3.189;3.939],2e-4,0);

results=struct('all48',all48,'table10_mesh',table10_mesh, ...
    'table10_smoothing',table10_smoothing);
fprintf('Reproduced Table 10 and all 48 Newton cases in %s\n',outdir)
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

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
