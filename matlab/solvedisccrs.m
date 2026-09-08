function [ux,uy,p,info]=solvedisccrs(tri,xnod,ynod,gam,dind,dval,opts)
% Stokes flow with the cavitation constraint p >= 0, discretized with the
% STABILIZED Crouzeix--Raviart element of Hansbo and Larson (M2AN 37 (2003)
% 63-72) and piecewise constant pressure.
%
% The plain CR element does not satisfy a discrete Korn inequality, so the
% strain based viscous form is unstable on it and one is forced to the full
% gradient form -- and hence, for the cavitation problem, to a constitutive
% law other than the deviatoric one.  Adding the (weakly consistent) jump
% penalty of that paper, eq. (26),
%
%     a(u,v) = sum_T ( sigma(u), eps(v) )_T
%              + 2*mu*gamma1 * sum_E ( h_E^{-1} [u], [v] )_E ,
%
% restores the discrete Korn inequality, and the deviatoric law may be used
% after all.  On the CR space the Nitsche consistency terms of the underlying
% discontinuous Galerkin method drop out identically: sigma(u) is elementwise
% constant for piecewise linears and the CR jump has vanishing mean on every
% edge, so ( <n.sigma(u)>, [v] )_E = 0.
%
% JUMP EVALUATION.  On edge j of element T (opposite vertex j) with global
% endpoints A and B,
%      u|_T = u_j + s*( u_alpha - u_beta ),   s = 1-2*lambda_alpha in [-1,1],
% where alpha, beta are the local indices of A, B, which are also the local
% indices of the other two edges.  The mean term u_j is the shared CR degree
% of freedom and cancels from the jump, leaving
%      [u] = s*d ,   d = (u_alpha-u_beta)|_{T+} - (u_alpha-u_beta)|_{T-} ,
% so that  int_E [u].[v] ds = (|E|/3) d_u . d_v ,  a rank one form per
% component in the four adjacent edge unknowns.  On a Dirichlet edge the
% constant and varying parts decouple (int_E s ds = 0), so penalizing the
% varying part alone is consistent with imposing the mean strongly.
%
%   in : tri,xnod,ynod  P1 mesh
%        gam            augmentation parameter (default 100)
%        dind,dval      prescribed velocity dofs (2*e-1, 2*e for edge e)
%        opts .gamma1   jump penalty parameter (default 1, per the paper)
%             .lamfac   first Lame parameter as a multiple of mu:
%                       -2/3 deviatoric (default), -1 full gradient form
%             .solver   'as' active set (default) or 'newton' on the smoothed
%                       operator phi_s(w) = w/2 + sqrt(w^2/4+s)
%             .s        smoothing parameter for 'newton' (default 1e-8)
%             .diredge  logical, one per edge: penalize this boundary edge
%                       (Dirichlet part).  Default: every boundary edge whose
%                       velocity dofs are all prescribed.
%             .dirichlet_components optional nedge-by-2 logical mask: boundary
%                       penalty on prescribed Cartesian components only.
%             .load     optional assembled velocity load vector (default zero).
%             .solver='qp' convex QP with checked working-set refinement.
%             .qp_tolerance .kkt_tolerance apply to the QP branch.
%             .reference_u .reference_p optional exact vectors for residual audit.
%             .verbose

if(nargin < 4 || isempty(gam)), gam=100; end
if(nargin < 7), opts=struct; end
if(~isfield(opts,'gamma1')),  opts.gamma1=1;      end
if(~isfield(opts,'lamfac')),  opts.lamfac=-2/3;   end
if(~isfield(opts,'solver')),  opts.solver='as';   end
if(~isfield(opts,'s')),       opts.s=1e-8;        end
if(~isfield(opts,'verbose')), opts.verbose=1;     end

nele=size(tri,1);my=1;
lam=opts.lamfac*my;
dmat=[lam+2*my, lam, 0; lam, lam+2*my, 0; 0, 0, my];

% ---------------- edges ----------------
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');
edg=reshape(ic,nele,3);
nedge=size(eu,1);
xe=(xnod(eu(:,1))+xnod(eu(:,2)))/2;
ye=(ynod(eu(:,1))+ynod(eu(:,2)))/2;
onbnd=(accumarray(ic,1)==1);
neqU=2*nedge;neq=neqU+nele;

% ---------------- elementwise forms ----------------
[ra,ca,va,ua]=assemble(nele*36);
rb=zeros(nele*6,1);cb=rb;vb=rb;ub=0;ar=zeros(nele,1);
for iel=1:nele
    iv=tri(iel,:);xc=xnod(iv);yc=ynod(iv);
    [~,fix,fiy,A]=basescalar(sum(xc)/3,sum(yc)/3,xc,yc);
    psix=-2*fix;psiy=-2*fiy;ar(iel)=A;
    Be=zeros(3,6);
    Be(1,1:2:end)=psix';
    Be(2,2:2:end)=psiy';
    Be(3,1:2:end)=psiy';Be(3,2:2:end)=psix';
    ael=A*(Be'*dmat*Be);
    bel=zeros(1,6);bel(1:2:end)=A*psix';bel(2:2:end)=A*psiy';
    ie=edg(iel,:);eqs=zeros(6,1);eqs(1:2:end)=2*ie-1;eqs(2:2:end)=2*ie;
    [ra,ca,va,ua]=assemble(ael,eqs,ra,ca,va,ua);
    rb(ub+1:ub+6)=iel;cb(ub+1:ub+6)=eqs;vb(ub+1:ub+6)=bel;ub=ub+6;
end

% ---------------- which boundary edges carry the penalty ----------------
if(nargin>=5 && ~isempty(dind))
    isD=false(nedge,1);
    pres=false(neqU,1);pres(dind(:))=true;
    for e=1:nedge
        if(onbnd(e) && pres(2*e-1) && pres(2*e)), isD(e)=true; end
    end
else
    isD=onbnd;
end
if(isfield(opts,'diredge') && ~isempty(opts.diredge)), isD=opts.diredge(:); end
dircomp=repmat(isD,1,2);
if(isfield(opts,'dirichlet_components'))
    dircomp=logical(opts.dirichlet_components);
    assert(isequal(size(dircomp),[nedge,2]));
    prescomp=reshape(pres,2,[])';
    assert(all(~dircomp(:)|prescomp(:)), ...
        'Boundary penalty requires prescribed component means');
end

% ---------------- jump penalty ----------------
% element areas adjacent to each edge, for h|_E = (|T+|+|T-|)/(2|E|)
elen=sqrt((xnod(eu(:,1))-xnod(eu(:,2))).^2+(ynod(eu(:,1))-ynod(eu(:,2))).^2);
arsum=accumarray(ic,repmat(ar,3,1),[nedge,1]);
hE=arsum./(2*elen);

% edge -> adjacent elements and their local edge index, built once
e2el=zeros(nedge,2);e2loc=zeros(nedge,2);cnt=zeros(nedge,1);
for j=1:3
    for iel=1:nele
        e=edg(iel,j);cnt(e)=cnt(e)+1;
        e2el(e,cnt(e))=iel;e2loc(e,cnt(e))=j;
    end
end

np=0;rp=zeros(nedge*2*16,1);cp=rp;vp=rp;
for e=1:nedge
    if(onbnd(e) && ~any(dircomp(e,:))), continue, end     % Neumann edge: no penalty
    A=eu(e,1);Bv=eu(e,2);
    dofs=[];coef=[];sgn=[1,-1];
    for side=1:cnt(e)
        iel=e2el(e,side);j=e2loc(e,side);iv=tri(iel,:);
        oth=setdiff(1:3,j);
        al=oth(iv(oth)==A);be=oth(iv(oth)==Bv);
        if(isempty(al)||isempty(be)), dofs=[];break, end
        dofs=[dofs,edg(iel,al),edg(iel,be)];
        coef=[coef,sgn(side),-sgn(side)];
    end
    if(isempty(dofs)), continue, end
    w=2*my*opts.gamma1*(elen(e)/3)/hE(e);
    blk=w*(coef'*coef);n=numel(dofs);
    for comp=0:1
        if(onbnd(e) && ~dircomp(e,comp+1)), continue, end
        eqs=2*dofs-1+comp;
        [I,J]=ndgrid(eqs,eqs);
        rp(np+1:np+n*n)=I(:);cp(np+1:np+n*n)=J(:);vp(np+1:np+n*n)=blk(:);
        np=np+n*n;
    end
end
K=sparse([ra(1:ua);rp(1:np)],[ca(1:ua);cp(1:np)],[va(1:ua);vp(1:np)],neqU,neqU);
B=sparse(rb(1:ub),cb(1:ub),vb(1:ub),nele,neqU);
M=spdiags(ar,0,nele,nele);Mi=spdiags(1./ar,0,nele,nele);

% ---------------- boundary conditions ----------------
if(nargin < 5 || isempty(dind))
    b1=find(onbnd & (abs(ye-max(ynod))<1e-9 | abs(ye-min(ynod))<1e-9));
    b2=find(onbnd & abs(xe-min(xnod))<1e-9);
    dind=[2*b1-1;2*b1;2*b2-1;2*b2];
    dval=[zeros(2*length(b1),1);ye(b2).*(1-ye(b2));zeros(length(b2),1)];
end
[ind,ia]=unique(dind(:),'first');bcval=dval(ia);
int=setdiff((1:neq)',ind);
ee=(1:nele)';

force=zeros(neqU,1);
if(isfield(opts,'load')),force=opts.load(:);assert(numel(force)==neqU);end
if(strcmp(opts.solver,'qp'))
  [u,p,act,qpinfo]=solve_cavitation_qp(K,B,ar,ind,bcval,force,opts);
  conv=qpinfo.converged;nit=qpinfo.iterations;
elseif(strcmp(opts.solver,'newton'))
  s=opts.s;
  u=zeros(neqU,1);u(ind)=bcval;p=zeros(nele,1);
  R=crresid(u,p,K,B,ar,gam,s,force);r0=max(norm(R(int)),1);nit=0;conv=false;
  for it=1:60
     nr=norm(R(int));
     if(nr<1e-10*r0), conv=true;nit=it-1;break, end
     w=gam*p-(B*u)./ar;[~,dph]=crsmax(w,s);
     D=spdiags(dph,0,nele,nele);
     J=[K+(1/gam)*(B'*D*Mi*B), -B'*D; -D*B, spdiags(gam*ar.*(dph-1),0,nele,nele)];
     del=zeros(neq,1);del(int)=J(int,int)\(-R(int));
     u=u+del(1:neqU);p=p+del(neqU+1:end);
     R=crresid(u,p,K,B,ar,gam,s,force);nit=it;
  end
  if(norm(R(int))<1e-10*r0), conv=true; end
  act=(gam*p-(B*u)./ar)>=0;
else
  u=zeros(neqU,1);p=zeros(nele,1);act=true(nele,1);actold=act;
  conv=false;nit=0;
  for it=1:100
    if(it==1), act=true(nele,1); else, act=(gam*p-(B*u)./ar>=0); end
    if(it>1 && isequal(act,actold)), conv=true;break, end
    actold=act;
    EA=sparse(ee(act),ee(act),1,nele,nele);EI=sparse(ee(~act),ee(~act),1,nele,nele);
    Smat=[K+(1/gam)*B'*EA*Mi*B, -B'*EA; -EA*B, -gam*EI*M];
    f=[force;zeros(nele,1)];fr=f-Smat(:,ind)*bcval;
    v=zeros(neq,1);v(ind)=bcval;v(int)=Smat(int,int)\fr(int);
    u=v(1:neqU);p=v(neqU+1:end);nit=it;
  end
end

ux=u(1:2:end);uy=u(2:2:end);D=(B*u)./ar;
checktol=1e-9*max([1,norm(p,inf),norm(D,inf)]);
signok=all(p>=-checktol) && all(D>=-checktol);
info=struct('iterations',nit,'converged',conv,'gamma',gam, ...
   'gamma1',opts.gamma1,'lamfac',opts.lamfac,'solver',opts.solver, ...
   'ncav',nnz(~act),'nele',nele,'cavfrac',nnz(~act)/nele,'act',act, ...
   'div',D,'edges',eu,'xe',xe,'ye',ye,'area',ar,'onbnd',onbnd, ...
   'compl',max(abs(p.*D)),'signok',signok,'checktol',checktol, ...
   'pmin',min(p),'pmax',max(p),'ndof',neq);
freeu=setdiff((1:neqU)',ind);
info.stationarity=norm(K(freeu,:)*u-force(freeu)-B(:,freeu)'*p,inf);
if(strcmp(opts.solver,'qp')),info.qp=qpinfo;end
if(isfield(opts,'reference_u') && isfield(opts,'reference_p'))
    info.reference_stationarity=norm(K(freeu,:)*opts.reference_u(:) ...
        -force(freeu)-B(:,freeu)'*opts.reference_p(:),inf);
    info.reference_divergence=norm(B*opts.reference_u(:)./ar,inf);
end
if(opts.verbose)
  fprintf(['solvedisccrs: %s, %d its (conv %d), gamma1 = %g, lamfac = %g\n' ...
     '   p in [%.4g, %.4g], min div u = %.3e, max|p*div u| = %.3e\n'], ...
     opts.solver,nit,conv,opts.gamma1,opts.lamfac,min(p),max(p),min(D),max(abs(p.*D)));
end
end

function [ph,dph]=crsmax(w,s)
r=sqrt(w.^2/4+s);ph=zeros(size(w));phm=ph;
q=w>=0;ph(q)=w(q)/2+r(q);ph(~q)=s./(r(~q)-w(~q)/2);
q=-w>=0;phm(q)=-w(q)/2+r(q);phm(~q)=s./(r(~q)+w(~q)/2);
dph=ph./(ph+phm);
end

function R=crresid(u,p,K,B,ar,gam,s,force)
w=gam*p-(B*u)./ar;[ph,~]=crsmax(w,s);lm=ph/gam;
R=[K*u-force-B'*lm; gam*ar.*(lm-p)];
end
