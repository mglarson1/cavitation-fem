function table05 = reproduce_stokes_3d_table(outdir)
% Reproduce manuscript Table 5: three-dimensional Stokes verification.
%
% The three coarser meshes use the direct saddle-point solver
% solvedisccrs3.  The 36x12x12 mesh uses solvedisccrs3_ip, which solves the
% same frozen systems by the iterated penalty method and needs about 6 GB of
% memory; the direct factorization does not fit in 16 GB.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

nx=[9;18;27;36];ny=nx/3;nz=ny;n=numel(nx);
elements=zeros(n,1);ndof=elements;iterations=elements;pmax=elements;
cavfrac=elements;max_abs_uz=elements;complementarity=elements;
solver=strings(n,1);
op=struct('gamma1',1,'lamfac',-2/3,'verbose',0);
for k=1:n
    [tet,x,y,z,dind,dval]=box3dcase(nx(k),ny(k),nz(k));
    if(k<n)
        [~,~,uz,p,info]=solvedisccrs3(tet,x,y,z,100,dind,dval,op);
        solver(k)="direct";
    else
        [~,~,uz,p,info]=solvedisccrs3_ip(tet,x,y,z,100,dind,dval,op);
        solver(k)="iterated penalty";
    end
    % P0 has no unique value on the mesh-aligned centreline.  The table
    % takes the maximum over the one-cell tube surrounding it.
    xc=mean(x(tet),2);yc=mean(y(tet),2);zc=mean(z(tet),2);
    tube=xc<2.5&abs(yc-0.5)<=1/ny(k)+1e-12& ...
        abs(zc-0.5)<=1/nz(k)+1e-12;
    pmax(k)=max(p(tube));
    elements(k)=size(tet,1);ndof(k)=info.ndof;
    iterations(k)=info.iterations;cavfrac(k)=info.cavfrac;
    max_abs_uz(k)=max(abs(uz));complementarity(k)=info.compl;
    assert(info.converged && info.signok)
end
table05=table(nx,ny,nz,elements,ndof,iterations,pmax,cavfrac, ...
    max_abs_uz,complementarity,solver);
writetable(table05,fullfile(outdir,'table09_stokes_3d_verification.csv'));

check_close('Table 5 elements',elements,[486;3888;13122;31104],0,0);
check_close('Table 5 degrees of freedom',ndof,[3780;28728;95256;223776],0,0);
check_close('Table 5 iterations',iterations,[2;3;3;3],0,0);
check_close('Table 5 peak pressures',pmax,[6.5538;6.0014;5.9806;5.9304],6e-5,0);
check_close('Table 5 cavity fractions',cavfrac,[0.0123;0.0183;0.0171;0.0186],6e-5,0);
check_close('Table 5 transverse velocities',max_abs_uz, ...
    [6.5e-3;2.9e-3;1.9e-3;1.4e-3],6e-5,0);
assert(max(complementarity(1:3))<5e-14 && complementarity(4)<5e-13)
fprintf('Reproduced Table 5 in %s\n',outdir)
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
