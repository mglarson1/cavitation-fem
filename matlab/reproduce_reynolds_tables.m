function results = reproduce_reynolds_tables(outdir)
% Reproduce manuscript Tables 1--4 and write machine-readable CSV files.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

delta=1;r=0.35;xp=1.5;yp=0.5;
rho=@(x,y) ((x-xp).^2+(y-yp).^2)/r^2;
dfun=@(x,y) 1+delta*exp(-rho(x,y));
dxfun=@(x,y) -delta*(2*(x-xp)/r^2).*exp(-rho(x,y));
dyfun=@(x,y) -delta*(2*(y-yp)/r^2).*exp(-rho(x,y));
ffun=@(x,y) -dxfun(x,y);

%% Table 1: nodal Reynolds refinement
refs=(2:6)';n=numel(refs);
nodes=zeros(n,1);iterations=nodes;pmax=nodes;cavfrac=nodes;pmin=nodes;
complementarity=nodes;signok=false(n,1);
for k=1:n
    nx=3*2^refs(k);ny=2^refs(k);
    [tri,x,y]=meshrect(3,1,nx,ny);
    evalc('[P,~,info]=solvereynolds(tri,x,y,dfun,ffun,100);');
    nodes(k)=numel(x);iterations(k)=info.iterations;pmax(k)=max(P);
    cavfrac(k)=info.cavfrac;pmin(k)=min(P);
    complementarity(k)=info.complementarity;signok(k)=info.signok;
    assert(info.converged && info.signok)
end
table01=table(refs,nodes,iterations,pmax,cavfrac,pmin, ...
    complementarity,signok);
writetable(table01,fullfile(outdir,'table01_reynolds_refinement.csv'));
check_close('Table 1 iteration counts',iterations,[3;4;8;13;22],0,0);
check_close('Table 1 peak pressures',pmax, ...
    [0.03999;0.04079;0.04181;0.04193;0.04194],6e-5,0);
check_close('Table 1 cavitated fractions',cavfrac, ...
    [0.231;0.311;0.353;0.383;0.398],6e-4,0);
assert(max(abs(pmin))<1e-20)

%% Table 2: manufactured obstacle solution
invh=[16;32;64;128;256];n=numel(invh);
l2=zeros(n,1);h1=zeros(n,1);iterations=zeros(n,1);
for k=1:n
    nx=2*invh(k);
    [tri,x,y]=meshrect(2,1,nx,1);x=x-1;
    dn=find(abs(x+1)<1e-12 | abs(x-1)<1e-12);
    done=@(x,y) ones(size(x));
    fman=@(x,y) (abs(x)<0.5).*(1-12*x.^2)+(abs(x)>=0.5).*(-2);
    evalc('[P,~,info]=solvereynolds(tri,x,y,done,fman,1,dn);');
    % The two traces carry equal and opposite diagonal-splitting errors.
    % Their mean is the one-dimensional P1 solution reported in the paper.
    [l2(k),h1(k)]=manufactured_line_errors(x,y,P);
    iterations(k)=info.iterations;
    assert(info.converged && info.signok)
end
rate_l2=[NaN;log(l2(1:end-1)./l2(2:end))/log(2)];
rate_h1=[NaN;log(h1(1:end-1)./h1(2:end))/log(2)];
table02=table(invh,l2,rate_l2,h1,rate_h1,iterations);
writetable(table02,fullfile(outdir,'table02_reynolds_convergence.csv'));
check_close('Table 2 L2 errors',l2, ...
    [3.15973508839871e-4;7.95505325765591e-5;1.99224101874956e-5; ...
     4.98277541281653e-6;1.24582964642452e-6],2e-13,0);
check_close('Table 2 H1 errors',h1, ...
    [1.60112682978433e-2;8.05295334145291e-3;4.03238766517915e-3; ...
     2.01693258581592e-3;1.00855863328477e-3],2e-13,0);
check_close('Table 2 iteration counts',iterations,[6;9;16;30;57],0,0);

%% Table 3: stability threshold
refs3=(2:4)';n=numel(refs3);
dof=zeros(n,1);lambda_max=zeros(n,1);gamma0_limit=zeros(n,1);
for k=1:n
    nx=3*2^refs3(k);ny=2^refs3(k);
    [tri,x,y]=meshrect(3,1,nx,ny);
    TR=triangulation(tri,x,y);
    bnod=double(abs(x)<1e-12|abs(x-3)<1e-12|abs(y)<1e-12|abs(y-1)<1e-12);
    [t6,x6,y6]=sqtriKb(tri,x,y,TR.edges,bnod);
    [K,H]=reynolds_stability_matrices(t6,x6,y6,dfun,dxfun,dyfun);
    E=[t6(:,[1,3]),t6(:,2);t6(:,[3,6]),t6(:,5);t6(:,[6,1]),t6(:,4)];
    [~,~,ic]=unique(sort(E(:,1:2),2),'rows');
    sel=ismember(ic,find(accumarray(ic,1)==1));
    dnodes=unique([reshape(E(sel,1:2),[],1);E(sel,3)]);
    free=setdiff((1:numel(x6))',dnodes);
    lambda_max(k)=eigs(H(free,free),K(free,free),1,'largestreal');
    dof(k)=numel(free);gamma0_limit(k)=1/lambda_max(k);
end
table03=table(refs3,dof,lambda_max,gamma0_limit);
writetable(table03,fullfile(outdir,'table03_reynolds_stability_threshold.csv'));
check_close('Table 3 degrees of freedom',dof,[161;705;2945],0,0);
check_close('Table 3 eigenvalues',lambda_max,[231;290;334],1.0,0);

%% Table 4: nodal and stabilized Reynolds methods
refs4=(3:5)';n=numel(refs4);
elements=zeros(n,1);nodal_iterations=elements;nodal_pmax=elements;
stabilized_iterations=elements;stabilized_pmax=elements;
stabilized_pmin=elements;difference=elements;
for k=1:n
    nx=3*2^refs4(k);ny=2^refs4(k);
    [tri,x,y]=meshrect(3,1,nx,ny);TR=triangulation(tri,x,y);
    bnod=double(abs(x)<1e-12|abs(x-3)<1e-12|abs(y)<1e-12|abs(y-1)<1e-12);
    [t6,x6,y6,b6]=sqtriKb(tri,x,y,TR.edges,bnod);
    evalc('[P1,~,i1]=solvereynolds(tri,x,y,dfun,ffun,100);');
    evalc('[P2,i2]=solvereynolds2(t6,x6,y6,dfun,dxfun,dyfun,1e-3);');
    elements(k)=size(tri,1);nodal_iterations(k)=i1.iterations;
    nodal_pmax(k)=max(P1);stabilized_iterations(k)=i2.iterations;
    stabilized_pmax(k)=max(P2);stabilized_pmin(k)=min(P2);
    difference(k)=norm(P2(1:numel(P1))-P1,inf);
    assert(i1.converged && i1.signok && i2.converged)
end
table04=table(refs4,elements,nodal_iterations,nodal_pmax, ...
    stabilized_iterations,stabilized_pmax,stabilized_pmin,difference);
writetable(table04,fullfile(outdir,'table04_reynolds_method_comparison.csv'));
check_close('Table 4 nodal iterations',nodal_iterations,[4;8;13],0,0);
check_close('Table 4 stabilized iterations',stabilized_iterations,[14;24;44],0,0);
check_close('Table 4 nodal peaks',nodal_pmax, ...
    [0.0407945;0.0418111;0.0419265],6e-7,0);
check_close('Table 4 stabilized peaks',stabilized_pmax, ...
    [0.0419143;0.0419537;0.0419507],6e-7,0);
check_close('Table 4 differences',difference, ...
    [6.1e-4;2.2e-4;3.9e-5],6e-6,0);

results=struct('table01',table01,'table02',table02, ...
    'table03',table03,'table04',table04);
fprintf('Reproduced Tables 1--4 in %s\n',outdir)
end

function [l2,h1]=manufactured_line_errors(x,y,P)
% Average the two traces of the one-element-wide strip, then integrate the
% resulting piecewise-linear one-dimensional pressure to high accuracy.
yval=unique(y);
assert(numel(yval)==2)
j0=find(abs(y-yval(1))<1e-12);j1=find(abs(y-yval(2))<1e-12);
[xs,o0]=sort(x(j0));[xs1,o1]=sort(x(j1));
assert(max(abs(xs-xs1))<1e-14)
ps=0.5*(P(j0(o0))+P(j1(o1)));
pex=@(z) (abs(z)<0.5).*(0.25-z.^2).^2;
dpex=@(z) (abs(z)<0.5).*(-4*z.*(0.25-z.^2));
l2sq=0;h1sq=0;
for j=1:numel(xs)-1
    a=xs(j);b=xs(j+1);slope=(ps(j+1)-ps(j))/(b-a);
    ph=@(z) ps(j)+slope*(z-a);
    l2sq=l2sq+integral(@(z) (pex(z)-ph(z)).^2,a,b, ...
        'AbsTol',1e-18,'RelTol',1e-13);
    h1sq=h1sq+integral(@(z) (dpex(z)-slope).^2,a,b, ...
        'AbsTol',1e-18,'RelTol',1e-13);
end
l2=sqrt(l2sq);h1=sqrt(h1sq);
end

function [K,H]=reynolds_stability_matrices(tri,x,y,dfun,dxfun,dyfun)
nno=numel(x);[rk,ck,vk,uk]=assemble(size(tri,1)*36);
[rh,ch,vh,uh]=assemble(size(tri,1)*36);
for iel=1:size(tri,1)
    iv=tri(iel,:);ivp=iv([1,3,6]);xc=x(ivp);yc=y(ivp);
    Kel=zeros(6);Hel=zeros(6);hT2=2*polyarea(xc,yc);
    [gx,gy,gw]=trigauc(xc,yc,5);
    for q=1:numel(gw)
        [~,fix,fiy,fixx,fiyy,~,A]=baseq2(gx(q),gy(q),xc,yc);
        d=dfun(gx(q),gy(q));dx=dxfun(gx(q),gy(q));dy=dyfun(gx(q),gy(q));
        Lv=d^3*(fixx+fiyy)+3*d^2*(dx*fix+dy*fiy);
        w=gw(q)*A;
        Kel=Kel+w*d^3*(fix*fix'+fiy*fiy');
        Hel=Hel+w*hT2*(Lv*Lv');
    end
    [rk,ck,vk,uk]=assemble(Kel,iv,rk,ck,vk,uk);
    [rh,ch,vh,uh]=assemble(Hel,iv,rh,ch,vh,uh);
end
K=sparse(rk(1:uk),ck(1:uk),vk(1:uk),nno,nno);
H=sparse(rh(1:uh),ch(1:uh),vh(1:uh),nno,nno);
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);
lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
