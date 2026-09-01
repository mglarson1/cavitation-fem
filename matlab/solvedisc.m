function [ux,uy,p,info]=solvedisc(nodes,xnod,ynod,nnop,bnod,gam,dind,dval,opts)
global gammazero
% Stokes flow with the cavitation constraint p >= 0, solved with the
% augmented Lagrangian formulation (aug1h)-(aug2h) of the paper and
% Taylor-Hood P2/P1 elements.
%
% The nonlinearity is [gamma*p - div u]_+ .  It is resolved by a
% primal-dual active set (semismooth Newton) fixed point iteration:
% given (u_k,p_k) the sign of  gamma*p - div u  is frozen at every
% Gauss point, the resulting linear system is solved for (u_{k+1},p_{k+1}),
% and the process is repeated until the active set stops changing.
% An unchanged active set reproduces the same linear system and is therefore
% an exact fixed point.  Finitely many possible sets do not by themselves
% exclude cycling, so the iteration also has a pressure tolerance and a
% maximum-iteration safeguard.
%
% Local numbering of the quadratic triangle:
%
%     6
%     | \
%     4   5
%     |    \
%     1--2--3
%
% Vertices (pressure nodes) are local nodes 1, 3, 6.
%
%   in : nodes  nele x 6 connectivity of the P2 mesh
%        xnod,ynod   coordinates of all P2 nodes
%        nnop        number of pressure (vertex) nodes
%        bnod        boundary flag
%        gam         (optional) augmentation parameter, default 100
%        dind        (optional) prescribed velocity degrees of freedom,
%                    numbered 2*node-1 for u_x and 2*node for u_y.  The
%                    default is the channel of maintri.m: no slip on the
%                    top and bottom walls and a parabolic inflow profile on
%                    the left, with the right hand end left traction free.
%        dval        (optional) the values prescribed at dind
%        opts.lamfac (optional) first Lame parameter divided by mu;
%                    default -2/3 for the deviatoric law, use 0 for the
%                    customary incompressible Stokes stress comparison
%
%  The viscous law is the deviatoric one,
%      sigma = 2*mu*(eps(u) - (1/3)*div(u)*I) - p*I ,
%  the Newtonian law with zero bulk viscosity, with the factor 1/3 of the
%  physical, three dimensional stress kept also in this plane computation.
%  Taking the three dimensional trace gives tr(sigma)/3 = -p, so p is the
%  mean normal pressure and p >= 0 has its intended mechanical meaning.
%  out : ux,uy       velocity, nno x 1 each
%        p           pressure, nnop x 1
%        info        struct with iteration diagnostics

if(nargin < 6 || isempty(gam)), gam = 100; end
if(nargin < 9), opts=struct; end
if(~isfield(opts,'lamfac')), opts.lamfac=-2/3; end
gammazero = gam;

maxite = 50;      % safeguard against cycling of the active set
ptol   = 1e-10;   % secondary, pressure based, stopping tolerance

bweps=zeros(3,12);
bw=zeros(2,12);
div=zeros(1,12);ieqs=div';

[nele,~]=size(nodes);
nno=length(xnod);

neq=2*nno+nnop;
my=1;
% The deviatoric law in Lame form.  The factor is 3 also in this plane
% (two dimensional) computation: the physical stress is three dimensional,
% and only lambda = -2*mu/3 makes p the mean of the three normal stresses.
% It also keeps the form coercive: with -2*mu/2 the two dimensional
% trace-free Korn inequality would be needed, and it does not hold.
lambda=opts.lamfac*my;
dmat=[lambda+2*my, lambda, 0; lambda, lambda+2*my, 0; 0, 0, my];

poldeg=3;

%-----------------------------------------------------------------------
%  Dirichlet data.  Independent of the iterate, so set up once.
%  Note the unique(): the two upstream corners belong both to the
%  wall sets and to the inflow set, and a repeated index would make
%  f = f - Smat(:,ind)*bcval subtract that column twice.
%-----------------------------------------------------------------------
if(nargin >= 8 && ~isempty(dind))
    ind=dind(:);bcval=dval(:);
else
    ind=[];bcval=[];

    ind1=find(ynod == max(ynod));                   % fixed top wall
    ind=[ind;2*ind1-1;2*ind1];
    bcval=[bcval;zeros(size(ind1));zeros(size(ind1))];

    ind1=find(ynod == min(ynod));                   % fixed bottom wall
    ind=[ind;2*ind1-1;2*ind1];
    bcval=[bcval;zeros(size(ind1));zeros(size(ind1))];

    ind1=find(xnod == min(xnod));                   % parabolic inflow
    ind=[ind;2*ind1-1;2*ind1];
    bcval=[bcval;ynod(ind1).*(1-ynod(ind1));zeros(size(ind1))];
end

[ind,ia]=unique(ind(:),'first');bcval=bcval(ia);
int=setdiff((1:neq)',ind);

%-----------------------------------------------------------------------
%  Active set iteration
%-----------------------------------------------------------------------
ux=zeros(nno,1);uy=ux;p=zeros(nnop,1);

[~,~,gv0]=trigauc(xnod(nodes(1,[1,3,6])),ynod(nodes(1,[1,3,6])),poldeg);
ngau=length(gv0);
act=false(nele,ngau);actold=act;

converged=false;pnorm=NaN;pcomp=1;nit=0;
for iite=1:maxite

    f=zeros(neq,1);
    nsize=nele*15*15;              % 15 = 12 velocity + 3 pressure dofs
    row=zeros(nsize,1);
    col=row;
    val=row;
    up=0;

    for iel=1:nele

        sele=zeros(12,12);fele=zeros(12,1);selep=zeros(3,12);selem=zeros(3);

        iv=nodes(iel,:);
        ivp=iv([1,3,6]);
        xc=xnod(ivp);
        yc=ynod(ivp);

        [gcx,gcy,gv]=trigauc(xc,yc,poldeg);
        for igau=1:length(gv)
            xm=gcx(igau);ym=gcy(igau);
            [fip,~,~,~]=basis(xm,ym,xc,yc);
            [fi,fix,fiy,detj]=baseq(xm,ym,xc,yc);

            bweps(1,1:2:end)=fix';
            bweps(2,2:2:end)=fiy';
            bweps(3,1:2:end)=fiy';
            bweps(3,2:2:end)=fix';
            bwsig=dmat*bweps;

            bw(1,1:2:end)=fi';
            bw(2,2:2:end)=fi';
            div(1:2:end)=fix';div(2:2:end)=fiy';

%             uxl=dot(fi,ux(iv));uyl=dot(fi,uy(iv));
%             bconv(1,1:2:end)=uxl*fix'+uyl*fiy';
%             bconv(2,2:2:end)=uxl*fix'+uyl*fiy';
%             sele=sele+gv(igau)*detj*bw'*bconv;

            fele=fele+0*gv(igau)*detj*bw'*floa(xm,ym);
            sele=sele+gv(igau)*detj*(bwsig)'*bweps;

            ploc=dot(fip,p(ivp));
            divloc=dot(fix,ux(iv))+dot(fiy,uy(iv));
            g=gammazero;
            test=(g*ploc-divloc);
            act(iel,igau)=(test >= 0);
            if(act(iel,igau))
                selep=selep-gv(igau)*detj*fip*div;
                sele=sele+gv(igau)*detj*(1/g)*(div)'*div;
            else
                selem=selem-gv(igau)*detj*g*(fip)*fip';
            end

        end

        ieqs(1:2:end)=2*iv-1;ieqs(2:2:end)=2*iv;
        ieqsp=ivp+2*nno;
        ivtot=[ieqs(:);ieqsp(:)];
        ele=[sele,selep';selep,selem];
        len=length(ivtot);
        X=ivtot(:,ones(1,len));
        Y=X';
        nn=len*len;
        lo=up+1;
        up=up+nn;
        row(lo:up)=X(:);
        col(lo:up)=Y(:);
        val(lo:up)=ele(:);

        f(ieqs)=f(ieqs)+fele;
    end

    % An unchanged active set reproduces the previous linear system, so
    % (ux,uy,p) from the last solve already is the fixed point.
    if(iite > 1 && isequal(act,actold))
        converged=true;
        break
    end
    actold=act;

    Smat=sparse(row(1:up),col(1:up),val(1:up),neq,neq);

    fr=f-Smat(:,ind)*bcval;
    fr=fr(int);
    u=full(Smat(int,int)\fr);
    v=zeros(neq,1);v(ind)=bcval;v(int)=u;
    ux=v(1:2:2*nno);uy=v(2:2:2*nno);

    nit=iite;
    pold=p;
    p=v((2*nno+1):end);
    if(iite == 1)
        pcomp=norm(p);
        if(pcomp == 0), pcomp=1; end
    end
    pnorm=norm(pold-p)/pcomp;
    fprintf('  iter %2d : active %6d/%6d gauss pts, |dp|/|p1| = %.3e\n', ...
        iite,nnz(act),numel(act),pnorm);
    if(pnorm < ptol)
        converged=true;
        break
    end
end

if(~converged)
    warning('solvedisc:noConvergence', ...
        'active set did not settle in %d iterations (last |dp|/|p1| = %.3e)', ...
        maxite,pnorm);
end

info=struct('iterations',nit,'converged',converged,'pnorm',pnorm, ...
    'gamma',gammazero,'nactive',nnz(act),'ngauss',numel(act), ...
    'cavelem',mean(~act,2), ...
    'pmin',min(p),'pmax',max(p),'lamfac',opts.lamfac,'ndof',neq);
fprintf(['solvedisc: %d iterations, gamma = %g, cavitated at %d of %d ' ...
    'gauss points, p in [%.4g, %.4g]\n'], ...
    nit,gammazero,numel(act)-nnz(act),numel(act),min(p),max(p));
