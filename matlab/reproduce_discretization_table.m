function table08 = reproduce_discretization_table(outdir)
% Reproduce manuscript Table 7: three Stokes discretizations.  Full-gradient
% CR on the two finer meshes is stopped by the sign rule (tol = 1e-12);
% there the active set changes only by elements with pressure and
% divergence at rounding level and never repeats.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

nx=[128;256;512];nz=[8;16;32];n=numel(nx);
methods={'proposed';'unstabilized';'taylor-hood'};nr=3*n;
mesh_nx=zeros(nr,1);mesh_nz=mesh_nx;method=cell(nr,1);ndof=mesh_nx;
iterations=mesh_nx;converged=false(nr,1);pmax=mesh_nx;cavity=mesh_nx;
pmin=mesh_nx;complementarity=nan(nr,1);
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
for k=1:n
    [tri,xh,zh,xm,zm,dind,dval,Hf]=pit_cr_mesh(nx(k),nz(k));

    [~,~,p,iP]=solvedisccrs(tri,xm,zm,100,dind,dval,opCR);
    mP=cr_metrics(tri,xm,zm,p,iP,Hf);
    j=3*(k-1)+1;
    [mesh_nx(j),mesh_nz(j),method{j}]=deal(nx(k),nz(k),methods{1});
    ndof(j)=iP.ndof;iterations(j)=iP.iterations;converged(j)=iP.converged;
    pmax(j)=mP.peak;cavity(j)=mP.cavity;pmin(j)=iP.pmin;
    complementarity(j)=iP.compl;

    tolU=[];if(k>1), tolU=1e-12; end
    evalc('[~,~,pU,iU]=solvedisccr(tri,xm,zm,100,dind,dval,tolU);');
    mU=cr_metrics(tri,xm,zm,pU,iU,Hf);
    j=j+1;
    [mesh_nx(j),mesh_nz(j),method{j}]=deal(nx(k),nz(k),methods{2});
    ndof(j)=iU.ndof;iterations(j)=iU.iterations;converged(j)=iU.converged;
    pmax(j)=mU.peak;cavity(j)=mU.cavity;pmin(j)=iU.pmin;
    complementarity(j)=iU.compl;

    mT=taylor_hood_metrics(nx(k),nz(k),Hf);
    j=j+1;
    [mesh_nx(j),mesh_nz(j),method{j}]=deal(nx(k),nz(k),methods{3});
    ndof(j)=mT.info.ndof;iterations(j)=mT.info.iterations;
    converged(j)=mT.info.converged;pmax(j)=mT.peak;cavity(j)=mT.cavity;
    pmin(j)=mT.info.pmin;
end
method=string(method);
table08=table(mesh_nx,mesh_nz,method,ndof,iterations,converged,pmax, ...
    cavity,pmin,complementarity);
writetable(table08,fullfile(outdir,'table08_discretization_comparison.csv'));

check_close('Table 7 degrees of freedom',ndof, ...
    [8464;8464;9899;33312;33312;38227;132160;132160;150179],0,0);
check_close('Table 7 proposed iterations',iterations(1:3:end),[8;9;10],0,0);
check_close('Table 7 Taylor--Hood iterations',iterations(3:3:end),[12;13;15],0,0);
assert(converged(2) && iterations(2)==9)
assert(all(converged([5,8])) && all(iterations([5,8])==9))
check_close('Table 7 peak pressures',pmax, ...
    [2.1537;1.9870;2.1556;2.1473;2.0934;2.1466;2.1451;2.1239;2.1452],6e-5,0);
check_close('Table 7 cavity lengths',cavity, ...
    [0.750;5.937;0.447;0.843;5.937;0.498;0.798;5.889;0.690],2e-4,0);
check_close('Table 7 Taylor--Hood undershoots',pmin(3:3:end), ...
    [-5.3e-4;-2.8e-4;-1.8e-4],6e-5,0);
assert(max(complementarity(1:3:end))<5e-13)
fprintf('Reproduced Table 7 in %s\n',outdir)
end

function [tri,xh,zh,xm,zm,dind,dval,Hf]=pit_cr_mesh(nx,nz)
L=24;xp=12;c=1;V=1;delta=1;rp=1;
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);xm=xh;zm=c-zh.*Hf(xh);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
sliding=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
shaped=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
end

function row=cr_metrics(tri,x,z,p,info,Hf)
L=24;c=1;xs=linspace(0,L,8001)';TR=triangulation(tri,x,z);
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
ps=nan(size(xs));ps(ok)=p(eid(ok));cav=false(size(xs));
cav(ok)=~info.act(eid(ok));w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(cav(w)));
end

function row=taylor_hood_metrics(nx,nz,Hf)
L=24;c=1;V=1;
[tri,xh,zh]=meshrect(L,1,nx,nz);TR=triangulation(tri,xh,zh);
[t6,x6,zh6,b6]=sqtriKb(tri,xh,zh,TR.edges,ones(size(xh)));
z6=c-zh6.*Hf(x6);nnop=numel(xh);
sliding=find(abs(zh6)<1e-9);shaped=find(abs(zh6-1)<1e-9);
dind=[2*sliding-1;2*sliding;2*shaped-1;2*shaped];
dval=[V*ones(numel(sliding),1);zeros(numel(sliding),1); ...
    zeros(numel(shaped),1);zeros(numel(shaped),1)];
opts=struct('lamfac',-2/3);
evalc('[~,~,p,iTH]=solvedisc(t6,x6,z6,nnop,b6,100,dind,dval,opts);');

xs=linspace(0,L,8001)';points=[xs,c-0.5*Hf(xs)];
TRm=triangulation(tri,xh,c-zh.*Hf(xh));eid=pointLocation(TRm,points);
ok=~isnan(eid);ps=nan(size(xs));
bc=cartesianToBarycentric(TRm,eid(ok),points(ok,:));
ps(ok)=sum(bc.*p(tri(eid(ok),:)),2);
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peak',max(ps(w)),'cavity',dx*nnz(ps(w)<=1e-8),'info',iTH);
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
