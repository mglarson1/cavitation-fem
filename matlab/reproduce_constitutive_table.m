function table07 = reproduce_constitutive_table(outdir)
% Reproduce manuscript Table 7: two constitutive laws on Taylor--Hood.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

cases=[0.5 8;1 4;1 2;1 1];n=size(cases,1);
delta=cases(:,1);radius=cases(:,2);steepness=delta./radius;
peak_deviatoric=zeros(n,1);peak_customary=peak_deviatoric;
cavity_deviatoric=peak_deviatoric;cavity_customary=peak_deviatoric;
iterations_deviatoric=peak_deviatoric;iterations_customary=peak_deviatoric;
for k=1:n
    A=taylor_hood_pit_metrics(delta(k),radius(k),256,16,-2/3);
    B=taylor_hood_pit_metrics(delta(k),radius(k),256,16,0);
    peak_deviatoric(k)=A.peak;peak_customary(k)=B.peak;
    cavity_deviatoric(k)=A.cavity;cavity_customary(k)=B.cavity;
    iterations_deviatoric(k)=A.info.iterations;
    iterations_customary(k)=B.info.iterations;
    assert(A.info.converged && B.info.converged)
end
peak_difference_pct=100*abs(peak_customary-peak_deviatoric)./peak_deviatoric;
cavity_ratio=cavity_customary./cavity_deviatoric;
table07=table(delta,radius,steepness,peak_deviatoric,peak_customary, ...
    peak_difference_pct,cavity_deviatoric,cavity_customary,cavity_ratio, ...
    iterations_deviatoric,iterations_customary);
writetable(table07,fullfile(outdir,'table07_constitutive_law_comparison.csv'));

check_close('Table 7 deviatoric peaks',peak_deviatoric, ...
    [3.3738;5.4539;3.7559;2.1466],6e-5,0);
check_close('Table 7 customary peaks',peak_customary, ...
    [3.3566;5.4336;3.7348;2.1253],6e-5,0);
check_close('Table 7 deviatoric cavities',cavity_deviatoric, ...
    [1.935;2.706;1.632;0.498],2e-4,0);
check_close('Table 7 customary cavities',cavity_customary, ...
    [2.106;2.847;2.073;1.113],2e-4,0);
fprintf('Reproduced Table 7 in %s\n',outdir)
end

function row=taylor_hood_pit_metrics(delta,rp,nx,nz,lamfac)
L=24;xp=L/2;c=1;V=1;gam=100;
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);TR=triangulation(tri,xh,zh);
[t6,x6,zh6,b6]=sqtriKb(tri,xh,zh,TR.edges,ones(size(xh)));
z6=c-zh6.*Hf(x6);nnop=numel(xh);
sliding=find(abs(zh6)<1e-9);shaped=find(abs(zh6-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
opts=struct('lamfac',lamfac);
evalc('[~,~,p,iTH]=solvedisc(t6,x6,z6,nnop,b6,gam,dind,dval,opts);');

xs=linspace(0,L,8001)';points=[xs,c-0.5*Hf(xs)];
TRm=triangulation(tri,xh,c-zh.*Hf(xh));
eid=pointLocation(TRm,points);ok=~isnan(eid);ps=nan(size(xs));
bc=cartesianToBarycentric(TRm,eid(ok),points(ok,:));
pv=p(tri(eid(ok),:));ps(ok)=sum(bc.*pv,2);
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(ps(w)<=1e-8), ...
    'pmin',min(ps(w)),'info',iTH);
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
