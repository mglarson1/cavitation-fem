function [ux,uy,p,info]=solvedisccrn(tri,xnod,ynod,gam,s,dind,dval,opts)
% Stokes flow with the cavitation constraint p >= 0, Crouzeix--Raviart
% velocity and piecewise constant pressure, with the positive part replaced
% by the smoothed operator of the differentiable Nitsche contact paper,
%
%     phi_s(w) = w/2 + sqrt(w^2/4 + s) ,      phi_0 = max(0,.)
%
% and solved by Newton's method instead of a primal-dual active set strategy.
%
% The method is stationarity of the regularized form of the functional
% (augmin) of the cavitation paper,
%
%   F_s(u,p) = 1/2 a(u,u) - L(u) + (1/gam) int Phi_s(gam p - div u) - gam/2 int p^2
%
% with Phi_s' = phi_s.  Writing lam = (1/gam) phi_s(w), w = gam*p - div u,
% the residual is
%
%   R_u = K u - B' lam ,        R_p = gam M (lam - p)
%
% and, with D = diag(phi_s'(w)) which lies strictly between 0 and 1,
%
%   J = [ K + (1/gam) B' D Mi B     -B' D        ]
%       [        -D B             gam M (D - I)  ]
%
% If K is positive definite on the free velocity degrees of freedom, as it is
% for the boundary conditions used by the current drivers, this Jacobian has
% a positive definite (1,1) block and a negative definite (2,2) block. It is
% then quasidefinite and nonsingular for every s > 0 and every iterate; no
% inf-sup condition or active set enters. Away from switching points, letting
% s -> 0 sends D to the indicator of the active set and recovers the frozen
% system of solvedisccr.m. At a zero switching argument the limit is 1/2.
%
% The smoothed complementarity is p_T (div u)_T = s/gam elementwise, the
% central path, rather than p_T (div u)_T = 0.
%
%   in : tri,xnod,ynod   P1 mesh
%        gam             augmentation parameter (default 100)
%        s               smoothing parameter (default 1e-6)
%        dind,dval       prescribed velocity dofs (2*e-1, 2*e for edge e)
%        opts            .maxit .tol .verbose .damp
%  out : ux,uy,p,info

if(nargin < 4 || isempty(gam)), gam = 100; end
if(nargin < 5 || isempty(s)),   s   = 1e-6; end
if(nargin < 8), opts=struct; end
if(~isfield(opts,'maxit')),   opts.maxit=60;    end
if(~isfield(opts,'tol')),     opts.tol=1e-10;   end
if(~isfield(opts,'verbose')), opts.verbose=1;   end
if(~isfield(opts,'damp')),    opts.damp=1;      end   % allow backtracking

nele=size(tri,1);my=1;

% ---------------- edges ----------------
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');
edg=reshape(ic,nele,3);
nedge=size(eu,1);
xe=(xnod(eu(:,1))+xnod(eu(:,2)))/2;
ye=(ynod(eu(:,1))+ynod(eu(:,2)))/2;
onbnd=(accumarray(ic,1)==1);
neqU=2*nedge;neq=neqU+nele;

% ---------------- K, B, M ----------------
[ra,ca,va,ua]=assemble(nele*36);
rb=zeros(nele*6,1);cb=rb;vb=rb;ub=0;ar=zeros(nele,1);
for iel=1:nele
    iv=tri(iel,:);xc=xnod(iv);yc=ynod(iv);
    [~,fix,fiy,A]=basescalar(sum(xc)/3,sum(yc)/3,xc,yc);
    psix=-2*fix;psiy=-2*fiy;ar(iel)=A;
    ael=zeros(6);bel=zeros(1,6);
    ael(1:2:end,1:2:end)=my*A*(psix*psix'+psiy*psiy');
    ael(2:2:end,2:2:end)=my*A*(psix*psix'+psiy*psiy');
    bel(1:2:end)=A*psix';bel(2:2:end)=A*psiy';
    ie=edg(iel,:);eqs=zeros(6,1);eqs(1:2:end)=2*ie-1;eqs(2:2:end)=2*ie;
    [ra,ca,va,ua]=assemble(ael,eqs,ra,ca,va,ua);
    rb(ub+1:ub+6)=iel;cb(ub+1:ub+6)=eqs;vb(ub+1:ub+6)=bel;ub=ub+6;
end
K=sparse(ra(1:ua),ca(1:ua),va(1:ua),neqU,neqU);
B=sparse(rb(1:ub),cb(1:ub),vb(1:ub),nele,neqU);
Mi=spdiags(1./ar,0,nele,nele);

% ---------------- boundary conditions ----------------
if(nargin < 6 || isempty(dind))
    b1=find(onbnd & (abs(ye-max(ynod))<1e-9 | abs(ye-min(ynod))<1e-9));
    b2=find(onbnd & abs(xe-min(xnod))<1e-9);
    dind=[2*b1-1;2*b1;2*b2-1;2*b2];
    dval=[zeros(2*length(b1),1);ye(b2).*(1-ye(b2));zeros(length(b2),1)];
end
[ind,ia]=unique(dind(:),'first');bcval=dval(ia);
int=setdiff((1:neq)',ind);

% ---------------- Newton ----------------
u=zeros(neqU,1);u(ind)=bcval;      % ind refers to velocity dofs only
p=zeros(nele,1);
rhist=zeros(opts.maxit,1);nhist=0;nit=0;converged=false;nback=0;

    function [ph,dph]=smax(w,s)
        r=sqrt(w.^2/4+s);
        ph=zeros(size(w));
        pos=w>=0;
        ph(pos)=w(pos)/2+r(pos);
        ph(~pos)=s./(r(~pos)-w(~pos)/2);      % conjugate form, no cancellation
        phm=zeros(size(w));                   % phi_s(-w), same trick
        neg=-w>=0;
        phm(neg)=-w(neg)/2+r(neg);
        phm(~neg)=s./(r(~neg)+w(~neg)/2);
        dph=ph./(ph+phm);                     % = 1/2 + (w/4)/r, stably
    end

    function R=resid(u,p)
        w=gam*p-(B*u)./ar;
        [ph,~]=smax(w,s);
        lam=ph/gam;
        R=[K*u-B'*lam; gam*ar.*(lam-p)];
    end

R=resid(u,p);
r0=norm(R(int));if(r0==0),r0=1;end
for it=1:opts.maxit
    nrm=norm(R(int));
    rhist(it)=nrm;
    nhist=it;
    if(opts.verbose)
        fprintf('   newton %2d : |R| = %.4e  (rel %.3e)\n',it-1,nrm,nrm/r0);
    end
    if(nrm < opts.tol*max(1,r0))
        converged=true;nit=it-1;break
    end
    w=gam*p-(B*u)./ar;
    [~,dph]=smax(w,s);
    D=spdiags(dph,0,nele,nele);
    J=[K+(1/gam)*(B'*D*Mi*B), -B'*D; -D*B, spdiags(gam*ar.*(dph-1),0,nele,nele)];
    del=zeros(neq,1);
    del(int)=J(int,int)\(-R(int));
    % full step, with backtracking only if the residual would increase
    t=1;
    for ls=1:20
        un=u+t*del(1:neqU);pn=p+t*del(neqU+1:end);
        Rn=resid(un,pn);
        if(norm(Rn(int)) < nrm || opts.damp==0), break, end
        t=t/2;nback=nback+1;
    end
    u=un;p=pn;R=Rn;nit=it;
end
if(~converged)
    nrm=norm(R(int));
    if(nrm < opts.tol*max(1,r0)), converged=true; end
end

ux=u(1:2:end);uy=u(2:2:end);
Dv=(B*u)./ar;
w=gam*p-Dv;[ph,~]=smax(w,s);lam=ph/gam;
cavthreshold=sqrt(s/gam);
info=struct('iterations',nit,'converged',converged,'gamma',gam,'s',s, ...
    'res',norm(R(int)),'relres',norm(R(int))/r0,'rhist',rhist(1:max(nhist,1)), ...
    'nbacktrack',nback,'cavthreshold',cavthreshold, ...
    'ncav',nnz(p<cavthreshold),'nele',nele, ...
    'cavfrac',nnz(p<cavthreshold)/nele,'div',Dv,'lam',lam, ...
    'edges',eu,'xe',xe,'ye',ye,'area',ar,'compl',max(abs(p.*Dv)), ...
    'central_path_error',max(abs(p.*Dv-s/gam)), ...
    'pmin',min(p),'pmax',max(p));
if(opts.verbose)
    fprintf(['solvedisccrn: %d Newton steps (conv %d, %d backtracks), gamma = %g, s = %g\n' ...
        '   p in [%.4g, %.4g], min div u = %.3e, max|p*div u| = %.3e (s/gam = %.3e)\n'], ...
        nit,converged,nback,gam,s,min(p),max(p),min(Dv),max(abs(p.*Dv)),s/gam);
end
end
