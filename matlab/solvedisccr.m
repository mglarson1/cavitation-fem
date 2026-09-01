function [ux,uy,p,info]=solvedisccr(tri,xnod,ynod,gam,dind,dval)
% Stokes flow with the cavitation constraint p >= 0, discretized with the
% nonconforming Crouzeix--Raviart element for the velocity and piecewise
% constants for the pressure.
%
% *** COMPARISON METHOD, NOT THE PROPOSED STABILIZED METHOD. ***
% This pair makes the complementarity conditions conform, but it forces the
% full gradient viscous form, since the discrete Korn inequality fails for
% nonconforming piecewise linears, and the paper uses the deviatoric law
% instead.  The two cannot be had together; see the remark on conformity.
%
% WHY THIS PAIR.  With Taylor--Hood the divergence of the discrete velocity
% is a discontinuous piecewise linear while the pressure is a continuous
% one, so the two spaces do not match and the complementarity conditions
% can only be imposed after testing against Q^h; the discrete pressure then
% need not be nonnegative.  Here
%
%      div V^h  =  Q^h  =  piecewise constants ,
%
% so p_T >= 0, (div u)_T >= 0 and p_T (div u)_T = 0 hold exactly, element
% by element, and both constraint sets are approximated from within.
%
% The price is that V^h is not contained in [H^1]^n, and that the discrete
% Korn inequality fails for nonconforming piecewise linears.  The viscous
% term is therefore taken in the form mu*grad u : grad v rather than
% 2*mu*eps(u):eps(v).  The two differ by mu*(div u, div v) and so only on
% the cavitated set; which of them is appropriate there is a constitutive
% question about the two phase mixture and not one this element can settle.
%
% Local numbering: degree of freedom i sits at the midpoint of the edge
% opposite vertex i, and the basis functions are psi_i = 1 - 2*lambda_i
% with lambda_i the barycentric coordinates, so that grad psi_i =
% -2 grad lambda_i.
%
%   in : tri          nele x 3 connectivity
%        xnod,ynod    vertex coordinates
%        gam          (optional) augmentation parameter, default 100
%        dind         (optional) prescribed velocity degrees of freedom,
%                     numbered 2*e-1 and 2*e for edge e.  Default: no slip
%                     on the horizontal walls and a parabolic inflow on the
%                     left, as in maintri.m.
%        dval         (optional) the values prescribed at dind
%  out : ux,uy        velocity at the edge midpoints
%        p            pressure, one value per element
%        info         struct with diagnostics, including the edge list and
%                     midpoint coordinates

if(nargin < 4 || isempty(gam)), gam = 100; end

nele=size(tri,1);
my=1;

%-----------------------------------------------------------------------
%  Edges.  Edge k of an element is the one opposite its vertex k.
%-----------------------------------------------------------------------
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');
edg=reshape(ic,nele,3);
nedge=size(eu,1);
xe=(xnod(eu(:,1))+xnod(eu(:,2)))/2;
ye=(ynod(eu(:,1))+ynod(eu(:,2)))/2;
onbnd=(accumarray(ic,1)==1);

neqU=2*nedge;
neq=neqU+nele;

%-----------------------------------------------------------------------
%  Stiffness, divergence and element areas.  None depend on the active
%  set, so all are built once.
%-----------------------------------------------------------------------
[ra,ca,va,ua]=assemble(nele*36);
rb=zeros(nele*6,1);cb=rb;vb=rb;ub=0;
ar=zeros(nele,1);
for iel=1:nele
    iv=tri(iel,:);xc=xnod(iv);yc=ynod(iv);
    [~,fix,fiy,A]=basescalar(sum(xc)/3,sum(yc)/3,xc,yc);
    psix=-2*fix;psiy=-2*fiy;               % constant on the element
    ar(iel)=A;
    ael=zeros(6);bel=zeros(1,6);
    ael(1:2:end,1:2:end)=my*A*(psix*psix'+psiy*psiy');
    ael(2:2:end,2:2:end)=my*A*(psix*psix'+psiy*psiy');
    bel(1:2:end)=A*psix';bel(2:2:end)=A*psiy';
    ie=edg(iel,:);
    eqs=zeros(6,1);eqs(1:2:end)=2*ie-1;eqs(2:2:end)=2*ie;
    [ra,ca,va,ua]=assemble(ael,eqs,ra,ca,va,ua);
    rb(ub+1:ub+6)=iel;cb(ub+1:ub+6)=eqs;vb(ub+1:ub+6)=bel;ub=ub+6;
end
K=sparse(ra(1:ua),ca(1:ua),va(1:ua),neqU,neqU);
B=sparse(rb(1:ub),cb(1:ub),vb(1:ub),nele,neqU);   % (Bu)_T = int_T div u
M=spdiags(ar,0,nele,nele);
Mi=spdiags(1./ar,0,nele,nele);

%-----------------------------------------------------------------------
%  Boundary conditions
%-----------------------------------------------------------------------
if(nargin < 5 || isempty(dind))
    b1=find(onbnd & (abs(ye-max(ynod))<1e-9 | abs(ye-min(ynod))<1e-9));
    b2=find(onbnd & abs(xe-min(xnod))<1e-9);
    dind=[2*b1-1;2*b1;2*b2-1;2*b2];
    dval=[zeros(2*length(b1),1);ye(b2).*(1-ye(b2));zeros(length(b2),1)];
end
[ind,ia]=unique(dind(:),'first');bcval=dval(ia);
int=setdiff((1:neq)',ind);

%-----------------------------------------------------------------------
%  Active set iteration.  Active = full film, inactive = cavitated.
%-----------------------------------------------------------------------
maxite=100;
u=zeros(neqU,1);p=zeros(nele,1);
act=true(nele,1);actold=act;
converged=false;nit=0;pnorm=NaN;pcomp=1;
cavhistory=zeros(maxite,1);
ee=(1:nele)';
for iite=1:maxite
    if(iite == 1)
        act=true(nele,1);                       % start from Stokes flow
    else
        act=(gam*p-(B*u)./ar >= 0);
    end
    cavhistory(iite)=nnz(~act);
    if(iite > 1 && isequal(act,actold))
        converged=true;
        break
    end
    actold=act;

    EA=sparse(ee(act),ee(act),1,nele,nele);
    EI=sparse(ee(~act),ee(~act),1,nele,nele);
    Smat=[K+(1/gam)*B'*EA*Mi*B, -B'*EA; -EA*B, -gam*EI*M];
    f=zeros(neq,1);
    fr=f-Smat(:,ind)*bcval;
    v=zeros(neq,1);v(ind)=bcval;
    v(int)=Smat(int,int)\fr(int);

    nit=iite;pold=p;
    u=v(1:neqU);p=v(neqU+1:end);
    if(iite==1), pcomp=norm(p); if(pcomp==0),pcomp=1;end; end
    pnorm=norm(pold-p)/pcomp;
    fprintf('  iter %2d : cavitated %5d/%5d elements, |dp|/|p1| = %.3e\n', ...
        iite,nnz(~act),nele,pnorm);
end
if(~converged)
    warning('solvedisccr:noConvergence','active set did not settle in %d iterations',maxite);
end

ux=u(1:2:end);uy=u(2:2:end);
D=(B*u)./ar;
info=struct('iterations',nit,'converged',converged,'pnorm',pnorm, ...
    'gamma',gam,'ncav',nnz(~act),'nele',nele,'cavfrac',nnz(~act)/nele, ...
    'act',act,'div',D,'edges',eu,'xe',xe,'ye',ye,'onbnd',onbnd,'area',ar, ...
    'pmin',min(p),'pmax',max(p),'compl',max(abs(p.*D)), ...
    'ndof',neq,'cavhistory',cavhistory(1:max(nit,1)));
fprintf(['solvedisccr: %d iterations, gamma = %g, cavitated in %d of %d ' ...
    'elements, p in [%.4g, %.4g], min div u = %.3e, max |p*div u| = %.3e\n'], ...
    nit,gam,nnz(~act),nele,min(p),max(p),min(D),max(abs(p.*D)));
