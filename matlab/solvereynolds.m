function [P,lam,info]=solvereynolds(tri,xnod,ynod,dfun,ffun,gam,dnodes)
% Reynolds lubrication problem with the cavitation constraint P >= 0,
% solved with the augmented Lagrangian formulation of Section 2 of the
% paper and continuous piecewise linear elements.
%
% Strong form (paper, Section 2):
%
%   -div(d^3 grad P) - lambda = f   in Omega
%                           P = 0   on dOmega
%     P >= 0,  lambda >= 0,  lambda*P = 0
%
% with  P = p c^2/(6 mu V),  d = H/c,  f = -dd/dx.  The Kuhn--Tucker
% conditions are replaced by the equivalent statement
%
%   lambda = (1/gamma) [gamma*lambda - P]_+ ,
%
% turning the variational inequality into the variational equality: find
% (P,lambda) in H^1_0(Omega) x L_2(Omega) such that
%
%   (d^3 grad P, grad q) - (1/gamma)([gamma*lambda-P]_+, q) = (f,q)
%   ( (1/gamma)[gamma*lambda-P]_+ - lambda, mu )            = 0
%
% for all q in H^1_0 and mu in L_2.
%
% NOTE ON SIGNS.  The minus multiplier term agrees with the strong form
% -div(d^3 grad P) - lambda = f.  A nonnegative lambda acts as the source
% that lifts P to the obstacle, exactly as in the obstacle problem.
%
% DISCRETIZATION.  P and lambda are both continuous piecewise linears on
% the same mesh.  The L_2 pairings in which lambda is tested are evaluated
% by nodal (lumped) quadrature, m_i = int phi_i dOmega.  Two reasons:
%
%  1. For piecewise linears the pointwise constraint is exactly the nodal
%     one, since P_h attains its extrema at the nodes:
%     P_h >= 0 in Omega  <=>  P_i >= 0 for every node i.  Nodal quadrature
%     is therefore not a simplification of the constraint but a faithful
%     rendering of it.
%  2. Consistent (Gauss) quadrature with equal order spaces for P and
%     lambda is inf--sup unstable and the frozen systems come out singular;
%     this was verified before settling on the present form.
%
% The nonlinearity is resolved as in solvedisc.m, by a primal--dual active
% set (semismooth Newton) iteration: the sign of gamma*lambda - P is frozen
% at every node, the resulting linear system is solved, and the process is
% repeated until the active set stops changing.  With
% A = {i interior : gamma*lambda_i - P_i >= 0} (the cavitated nodes) and
% I its complement, the frozen system is
%
%   [ K + Da/gamma   -Da   ] [ P   ]   [ F ]
%   [   -Da        -gamma Di] [ lam ] = [ 0 ]
%
% with K the d^3--weighted stiffness matrix, Da = diag(m_i, i in A) and
% Di = diag(m_i, i in I).  It is symmetric, has the same block structure as
% the Stokes solver -- P here plays the role of div u and lambda that of p
% -- and reduces to the classical primal--dual active set method: P_i = 0
% on A, lambda_i = 0 on I, and the two are recovered from each other's
% equation.  Dirichlet nodes are kept out of A, where P = 0 is already
% imposed and no multiplier is needed.
%
% Because complementarity holds exactly at the nodes, gamma cancels from
% the method altogether.  On A, P_i = 0 and active-set membership gives
% lambda_i >= 0.  On I, lambda_i = 0 and inactive-set membership gives
% -P_i < 0, hence P_i > 0.  Ties are assigned to A.  The iteration is
% therefore the classical primal--dual active-set method, and neither the
% converged solution nor the iteration count depends on gamma.  Gamma is
% retained only to keep the correspondence with the formulation visible.
%
%   in : tri          nele x 3 connectivity (P1 triangles)
%        xnod,ynod    node coordinates
%        dfun         handle, d(x,y), the scaled film thickness
%        ffun         handle, f(x,y) = -dd/dx
%        gam          (optional) augmentation parameter, default 100
%        dnodes       (optional) nodes where P = 0, default free boundary
%  out : P            scaled pressure, nno x 1
%        lam          multiplier, nno x 1
%        info         struct with iteration diagnostics

if(nargin < 6 || isempty(gam)), gam = 100; end

nele=size(tri,1);
nno=length(xnod);
neq=2*nno;                        % P in 1:nno, lambda in nno+1:2*nno
poldeg=3;

if(nargin < 7 || isempty(dnodes))
    TR=triangulation(tri,xnod(:),ynod(:));
    dnodes=unique(reshape(freeBoundary(TR),[],1));
end
dnodes=dnodes(:);
isdir=false(nno,1);isdir(dnodes)=true;

%-----------------------------------------------------------------------
%  Stiffness, load and lumped mass.  None of these depend on the active
%  set, so they are built once.
%-----------------------------------------------------------------------
[row,col,val,up]=assemble(nele*9);
F=zeros(nno,1);m=zeros(nno,1);
for iel=1:nele
    iv=tri(iel,:);
    xc=xnod(iv);yc=ynod(iv);
    sele=zeros(3);fele=zeros(3,1);mele=zeros(3,1);
    [gcx,gcy,gv]=trigauc(xc,yc,poldeg);
    for igau=1:length(gv)
        xm=gcx(igau);ym=gcy(igau);
        [fi,fix,fiy,area]=basescalar(xm,ym,xc,yc);
        w=gv(igau)*area;
        sele=sele+w*dfun(xm,ym)^3*(fix*fix'+fiy*fiy');
        fele=fele+w*ffun(xm,ym)*fi;
        mele=mele+w*fi;
    end
    [row,col,val,up]=assemble(sele,iv,row,col,val,up);
    F(iv)=F(iv)+fele;
    m(iv)=m(iv)+mele;
end
K=sparse(row(1:up),col(1:up),val(1:up),nno,nno);

%-----------------------------------------------------------------------
%  Active set iteration
%-----------------------------------------------------------------------
maxite=400;   % the set recedes one node layer per pass, so long thin
              % domains need many more passes than compact ones
P=zeros(nno,1);lam=zeros(nno,1);
act=false(nno,1);actold=act;
converged=false;pnorm=NaN;pcomp=1;nit=0;

nn=(1:nno)';
for iite=1:maxite

    % The first pass is the unconstrained Reynolds problem, the standard
    % warm start for an active set iteration.  A cold start from
    % P = lam = 0 would instead declare the whole domain cavitated.
    if(iite == 1)
        act=false(nno,1);
    else
        act=(gam*lam-P >= 0) & ~isdir;
    end

    if(iite > 1 && isequal(act,actold))
        converged=true;
        break
    end
    actold=act;

    Da=sparse(nn(act),nn(act),m(act),nno,nno);
    Di=sparse(nn(~act),nn(~act),m(~act),nno,nno);
    Smat=[K+Da/gam, -Da; -Da, -gam*Di];

    f=[F;zeros(nno,1)];
    ind=dnodes;                            % P = 0 on the boundary
    bcval=zeros(size(ind));
    int=setdiff((1:neq)',ind);

    fr=f-Smat(:,ind)*bcval;
    fr=fr(int);
    v=zeros(neq,1);v(ind)=bcval;
    v(int)=full(Smat(int,int)\fr);

    nit=iite;
    Pold=P;
    P=v(1:nno);lam=v(nno+1:end);
    if(iite == 1)
        pcomp=norm(P);
        if(pcomp == 0), pcomp=1; end
    end
    pnorm=norm(Pold-P)/pcomp;
    fprintf('  iter %2d : cavitated %6d/%6d nodes, |dP|/|P1| = %.3e\n', ...
        iite,nnz(act),nno,pnorm);
end

if(~converged)
    warning('solvereynolds:noConvergence', ...
        'active set did not settle in %d iterations',maxite);
end

checktol=1e-9*max([1,norm(P,inf),norm(lam,inf)]);
signok=all(P>=-checktol) && all(lam>=-checktol);
info=struct('iterations',nit,'converged',converged,'pnorm',pnorm, ...
    'gamma',gam,'ncav',nnz(act),'nnodes',nno,'cavfrac',nnz(act)/nno, ...
    'act',act,'Pmin',min(P),'Pmax',max(P), ...
    'lammin',min(lam),'lammax',max(lam),'complementarity',max(abs(P.*lam)), ...
    'signok',signok,'checktol',checktol,'K',K,'F',F,'m',m,'isdir',isdir);
fprintf(['solvereynolds: %d iterations, gamma = %g, cavitated at %d of %d ' ...
    'nodes, P in [%.4g, %.4g], lambda in [%.4g, %.4g]\n'], ...
    nit,gam,nnz(act),nno,min(P),max(P),min(lam),max(lam));
