function [P,info]=solvereynolds2(tri,xnod,ynod,dfun,dxfun,dyfun,gam0,dnodes)
% Reynolds lubrication with cavitation, P >= 0, discretized with the
% multiplier--free stabilized method of the paper -- the functional
% (perturbed) and its Euler--Lagrange equations (FEM), (stab_form0) --
% using continuous piecewise quadratics.
%
% The multiplier is replaced elementwise by its strong residual,
% lambda = -div(d^3 grad P_h) - f, which turns the augmented Lagrangian
% into a minimisation over V_h alone.  Writing
%
%   L v      := div(d^3 grad v) = d^3 lap v + 3 d^2 grad d . grad v
%   Q_gam(v) := -v - gam * L v
%
% the method is: find P_h in V_h such that
%
%   a(P_h,v) + < [Q(P_h) - gam f]_+ / gam , Q(v) >_h
%            - gam < L P_h + f , L v >_h  =  (f,v)      for all v in V_h,
%
% with a(P,v) = (d^3 grad P, grad v) and <.,.>_h the elementwise L2 form.
%
% WHY P2.  For piecewise linears on straight triangles lap v_h = 0, so
% L v_h collapses to 3 d^2 grad d . grad v_h and the method loses the
% entire second derivative of the multiplier approximation.  The paper
% allows P_k for k >= 1, but the substitution only carries its intended
% meaning from k = 2 on.
%
% CHOICE OF gamma.  Freezing the sign of R = Q(P_h) - gam f gives the
% symmetric frozen system
%
%   a(P,v) - gam <L P, L v>_h + chi_A <Q(P),Q(v)>_h / gam  =  ...
%
% in which the middle term enters with a NEGATIVE sign.  Coercivity
% therefore requires it to be dominated by a(.,.).  The inverse estimate
% ||L v_h||_T <= C h_T^{-1} ||d^3 grad v_h||_T forces
%
%   gam_T = gam0 * h_T^2 ,   h_T = sqrt(2|T|),
%
% with gam0 small enough; a fixed, mesh independent gamma makes the system
% indefinite on any sufficiently fine mesh.  The paper calls gamma
% "arbitrary", which is true of the continuous formulation but not of the
% discrete method (perturbed).
%
% The admissible range was measured for this element by computing the
% largest generalized eigenvalue of <L.,L.>_h weighted by h_T^2 against
% a(.,.).  It is essentially mesh independent -- 231, 290, 334 on three
% successive refinements -- so the frozen systems are positive definite
% for
%
%   gam0 < 1/334 = 0.003   (approximately),
%
% which is the P2 inverse inequality constant C = 18 in disguise.  The
% default below sits a factor three inside that bound.  Above it the
% method diverges: the active set never settles and the pressure picks up
% large negative excursions.
%
% Unlike the nodal method in solvereynolds.m, this one is a consistent
% stabilized scheme: P_h >= 0 is imposed weakly, so small undershoots of
% order gam0 remain and vanish under refinement.
%
%   in : tri          nele x 6 connectivity (P2, local order 1..6 with
%                     vertices at 1,3,6 as produced by sqtriKb)
%        xnod,ynod    coordinates of all P2 nodes
%        dfun         handle, d(x,y), the scaled film thickness
%        dxfun,dyfun  handles, dd/dx and dd/dy.  The load is f = -dd/dx.
%        gam0         (optional) stabilization constant, default 1e-3
%        dnodes       (optional) nodes where P = 0, default the boundary
%  out : P            scaled pressure at the P2 nodes
%        info         struct with iteration diagnostics

if(nargin < 7 || isempty(gam0)), gam0 = 1e-3; end

nele=size(tri,1);
nno=length(xnod);
poldeg=5;

%-----------------------------------------------------------------------
%  Dirichlet nodes.  An edge of the P2 mesh is on the boundary when it
%  belongs to exactly one element; local edges are (1,3)-2, (3,6)-5 and
%  (6,1)-4 in the numbering produced by sqtriKb.
%-----------------------------------------------------------------------
if(nargin < 8 || isempty(dnodes))
    E=[tri(:,[1,3]),tri(:,2);tri(:,[3,6]),tri(:,5);tri(:,[6,1]),tri(:,4)];
    [~,~,ic]=unique(sort(E(:,1:2),2),'rows');
    sel=ismember(ic,find(accumarray(ic,1)==1));
    dnodes=unique([reshape(E(sel,1:2),[],1);E(sel,3)]);
end
ind=dnodes(:);
bcval=zeros(size(ind));
int=setdiff((1:nno)',ind);

maxite=50;
ptol=1e-10;
P=zeros(nno,1);
[~,~,gv0]=trigauc(xnod(tri(1,[1,3,6])),ynod(tri(1,[1,3,6])),poldeg);
ngau=length(gv0);
act=false(nele,ngau);actold=act;
converged=false;pnorm=NaN;pcomp=1;nit=0;

for iite=1:maxite

    [row,col,val,up]=assemble(nele*36);
    F=zeros(nno,1);

    for iel=1:nele
        iv=tri(iel,:);
        ivp=iv([1,3,6]);
        xc=xnod(ivp);yc=ynod(ivp);
        sele=zeros(6);fele=zeros(6,1);
        hT2=2*polyarea(xc,yc);              % h_T^2
        gam=gam0*hT2;

        [gcx,gcy,gv]=trigauc(xc,yc,poldeg);
        for igau=1:length(gv)
            xm=gcx(igau);ym=gcy(igau);
            [fi,fix,fiy,fixx,fiyy,~,area]=baseq2(xm,ym,xc,yc);
            w=gv(igau)*area;
            dl=dfun(xm,ym);dxl=dxfun(xm,ym);dyl=dyfun(xm,ym);
            fl=-dxl;                        % f = -dd/dx

            Lv=dl^3*(fixx+fiyy)+3*dl^2*(dxl*fix+dyl*fiy);
            Qv=-fi-gam*Lv;

            sele=sele+w*(dl^3*(fix*fix'+fiy*fiy')-gam*(Lv*Lv'));
            fele=fele+w*(fl*fi+gam*fl*Lv);

            % The first pass drops the [.]_+ term, which is the
            % unconstrained Reynolds problem with the same stabilization.
            if(iite == 1)
                act(iel,igau)=false;
            else
                act(iel,igau)=((dot(Qv,P(iv))-gam*fl) >= 0);
            end
            if(act(iel,igau))
                sele=sele+w*(1/gam)*(Qv*Qv');
                fele=fele+w*fl*Qv;
            end
        end

        [row,col,val,up]=assemble(sele,iv,row,col,val,up);
        F(iv)=F(iv)+fele;
    end

    nchg=nnz(act ~= actold);
    if(iite > 1 && nchg == 0)
        converged=true;                     % exact fixed point
        break
    end
    if(iite > 2 && pnorm < ptol)
        % The set still flips at a handful of Gauss points on the free
        % boundary, but the pressure has stopped moving.  Unlike the nodal
        % method, where P = 0 holds exactly on the active set, here the
        % test involves lap P_h and chatters harmlessly at points where
        % P_h is within round-off of the free boundary.
        converged=true;
        break
    end
    actold=act;

    Smat=sparse(row(1:up),col(1:up),val(1:up),nno,nno);
    fr=F-Smat(:,ind)*bcval;
    fr=fr(int);
    Pold=P;
    Pnew=zeros(nno,1);Pnew(ind)=bcval;
    Pnew(int)=full(Smat(int,int)\fr);
    P=Pnew;

    nit=iite;
    if(iite == 1)
        pcomp=norm(P);
        if(pcomp == 0), pcomp=1; end
        pnorm=1;
    else
        pnorm=norm(Pold-P)/pcomp;
    end
    fprintf('  iter %2d : cavitated %6d/%6d gauss pts, %4d changed, |dP|/|P1| = %.3e\n', ...
        iite,nnz(act),numel(act),nchg,pnorm);
end

if(~converged)
    warning('solvereynolds2:noConvergence', ...
        'active set did not settle in %d iterations',maxite);
end

info=struct('iterations',nit,'converged',converged,'pnorm',pnorm, ...
    'nchanged',nnz(act ~= actold),'gamma0',gam0,'ncav',nnz(act),'ngauss',numel(act), ...
    'cavfrac',nnz(act)/numel(act),'cavelem',mean(act,2), ...
    'Pmin',min(P),'Pmax',max(P));
fprintf(['solvereynolds2: %d iterations, gamma0 = %g, cavitated at %d of %d ' ...
    'gauss points, P in [%.4g, %.4g]\n'], ...
    nit,gam0,nnz(act),numel(act),min(P),max(P));
