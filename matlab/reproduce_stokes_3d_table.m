function table09 = reproduce_stokes_3d_table(outdir)
% Reproduce manuscript Table 9: three-dimensional Stokes verification.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

nx=[9;18;27];ny=[3;6;9];nz=ny;n=numel(nx);
elements=zeros(n,1);ndof=elements;iterations=elements;pmax=elements;
cavfrac=elements;max_abs_uz=elements;complementarity=elements;
op=struct('gamma1',1,'lamfac',-2/3,'verbose',0);
for k=1:n
    [tet,x,y,z]=meshbox(3,1,1,nx(k),ny(k),nz(k));
    [dind,dval]=channel3d_boundary(tet,x,y,z);
    [~,~,uz,p,info]=solvedisccrs3(tet,x,y,z,100,dind,dval,op);
    % P0 has no unique value on the mesh-aligned centerline.  The original
    % table takes the maximum over the one-cell tube surrounding it.
    xc=mean(x(tet),2);yc=mean(y(tet),2);zc=mean(z(tet),2);
    tube=xc<2.5&abs(yc-0.5)<=1/ny(k)+1e-12& ...
        abs(zc-0.5)<=1/nz(k)+1e-12;
    pmax(k)=max(p(tube));
    elements(k)=size(tet,1);ndof(k)=info.ndof;
    iterations(k)=info.iterations;cavfrac(k)=info.cavfrac;
    max_abs_uz(k)=max(abs(uz));complementarity(k)=info.compl;
    assert(info.converged && info.signok)
end
table09=table(nx,ny,nz,elements,ndof,iterations,pmax,cavfrac, ...
    max_abs_uz,complementarity);
writetable(table09,fullfile(outdir,'table09_stokes_3d_verification.csv'));

check_close('Table 9 elements',elements,[486;3888;13122],0,0);
check_close('Table 9 degrees of freedom',ndof,[3780;28728;95256],0,0);
check_close('Table 9 iterations',iterations,[2;3;3],0,0);
check_close('Table 9 peak pressures',pmax,[6.5538;6.0014;5.9806],6e-5,0);
check_close('Table 9 cavity fractions',cavfrac,[0.0123;0.0183;0.0171],6e-5,0);
check_close('Table 9 transverse velocities',max_abs_uz, ...
    [6.5e-3;2.9e-3;1.9e-3],6e-5,0);
assert(max(complementarity)<5e-13)
fprintf('Reproduced Table 9 in %s\n',outdir)
end

function [dind,dval]=channel3d_boundary(tet,x,y,z)
nele=size(tet,1);
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
assert(nele*4==numel(ic))
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
