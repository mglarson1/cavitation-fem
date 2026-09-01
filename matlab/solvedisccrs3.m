function [ux,uy,uz,p,info]=solvedisccrs3(tet,xnod,ynod,znod,gam,dind,dval,opts)
% Three dimensional Stokes flow with the cavitation constraint p >= 0,
% discretized with the STABILIZED Crouzeix-Raviart element of Hansbo and
% Larson (M2AN 37 (2003) 63-72) and piecewise constant pressure.  The direct
% counterpart of solvedisccrs.m, whose documentation applies verbatim with
% edges replaced by faces.
%
% On a tetrahedron the nonconforming linear basis is psi_i = 1 - 3*lambda_i,
% one degree of freedom per face, sitting at the face centroid, and equal to
% the mean of the trace over that face.  On face j (opposite vertex j) the
% trace is
%       u|_F = u_j + sum_{X in F} u_{alpha_X} phi_X ,   phi_X = 1 - 3 lambda_X ,
% where alpha_X is the local index of the face vertex X, and phi_X has
% vanishing mean over F.  The shared face degree of freedom u_j cancels from
% the jump, leaving
%       [u] = sum_X d_X phi_X ,   d_X = u_{alpha_X}|_{T+} - u_{alpha_X}|_{T-} ,
% and, since int_F phi_X phi_Y = |F| (3*delta_XY - 1)/4,
%       int_F [u].[v] = (|F|/4) d_u' (3 I - J) d_v ,
% with J the 3 by 3 matrix of ones.  The face parameter is
% h|_F = (|T+|+|T-|)/(2|F|), the volume-over-area analogue of the two
% dimensional definition.
%
%   in : tet         nele x 4 connectivity
%        xnod,..     vertex coordinates
%        gam         augmentation parameter (default 100; the method does not
%                    depend on it, which is the point of this pair)
%        dind,dval   prescribed velocity dofs, numbered 3*f-2, 3*f-1, 3*f
%                    for face f
%        opts .gamma1 jump penalty (default 1)
%             .lamfac first Lame parameter over mu (default -2/3)
%             .verbose

if(nargin < 5 || isempty(gam)), gam=100; end
if(nargin < 8), opts=struct; end
if(~isfield(opts,'gamma1')),  opts.gamma1=1;    end
if(~isfield(opts,'lamfac')),  opts.lamfac=-2/3; end
if(~isfield(opts,'verbose')), opts.verbose=1;   end

nele=size(tet,1);my=1;
lm=opts.lamfac*my;
dmat=zeros(6);dmat(1:3,1:3)=lm;dmat(1:3,1:3)=dmat(1:3,1:3)+2*my*eye(3);
dmat(4,4)=my;dmat(5,5)=my;dmat(6,6)=my;

% ---------------- faces: face j is opposite vertex j ----------------
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');
fac=reshape(ic,nele,4);
nface=size(fu,1);
onbnd=(accumarray(ic,1)==1);
neqU=3*nface;neq=neqU+nele;

xf=(xnod(fu(:,1))+xnod(fu(:,2))+xnod(fu(:,3)))/3;
yf=(ynod(fu(:,1))+ynod(fu(:,2))+ynod(fu(:,3)))/3;
zf=(znod(fu(:,1))+znod(fu(:,2))+znod(fu(:,3)))/3;

% face areas
v1=[xnod(fu(:,2))-xnod(fu(:,1)),ynod(fu(:,2))-ynod(fu(:,1)),znod(fu(:,2))-znod(fu(:,1))];
v2=[xnod(fu(:,3))-xnod(fu(:,1)),ynod(fu(:,3))-ynod(fu(:,1)),znod(fu(:,3))-znod(fu(:,1))];
cr=cross(v1,v2,2);
farea=0.5*sqrt(sum(cr.^2,2));

% ---------------- elementwise forms ----------------
[ra,ca,va,ua]=assemble(nele*144);
rb=zeros(nele*12,1);cb=rb;vb=rb;ub=0;vol=zeros(nele,1);
for iel=1:nele
    iv=tet(iel,:);xc=xnod(iv);yc=ynod(iv);zc=znod(iv);
    [~,lx,ly,lz,V]=basis3(sum(xc)/4,sum(yc)/4,sum(zc)/4,xc,yc,zc);
    px=-3*lx;py=-3*ly;pz=-3*lz;          % grad psi_i, constant on the element
    vol(iel)=V;
    Be=zeros(6,12);
    Be(1,1:3:end)=px';                    % eps_xx
    Be(2,2:3:end)=py';                    % eps_yy
    Be(3,3:3:end)=pz';                    % eps_zz
    Be(4,1:3:end)=py';Be(4,2:3:end)=px';  % gamma_xy
    Be(5,1:3:end)=pz';Be(5,3:3:end)=px';  % gamma_xz
    Be(6,2:3:end)=pz';Be(6,3:3:end)=py';  % gamma_yz
    ael=V*(Be'*dmat*Be);
    bel=zeros(1,12);
    bel(1:3:end)=V*px';bel(2:3:end)=V*py';bel(3:3:end)=V*pz';
    jf=fac(iel,:);eqs=zeros(12,1);
    eqs(1:3:end)=3*jf-2;eqs(2:3:end)=3*jf-1;eqs(3:3:end)=3*jf;
    [ra,ca,va,ua]=assemble(ael,eqs,ra,ca,va,ua);
    rb(ub+1:ub+12)=iel;cb(ub+1:ub+12)=eqs;vb(ub+1:ub+12)=bel;ub=ub+12;
end

% ---------------- which boundary faces carry the penalty ----------------
if(nargin>=6 && ~isempty(dind))
    pres=false(neqU,1);pres(dind(:))=true;
    isD=onbnd & pres(3*(1:nface)'-2) & pres(3*(1:nface)'-1) & pres(3*(1:nface)');
else
    isD=onbnd;
end

% ---------------- jump penalty ----------------
volsum=accumarray(ic,repmat(vol,4,1),[nface,1]);
hF=volsum./(2*farea);

f2el=zeros(nface,2);f2loc=zeros(nface,2);cnt=zeros(nface,1);
for j=1:4
    for iel=1:nele
        f=fac(iel,j);cnt(f)=cnt(f)+1;
        f2el(f,cnt(f))=iel;f2loc(f,cnt(f))=j;
    end
end

MF0=(3*eye(3)-ones(3))/4;                 % int_F phi_X phi_Y = |F| * MF0
np=0;maxn=nface*3*36*36;
rp=zeros(maxn,1);cp=rp;vp=rp;
for f=1:nface
    if(onbnd(f) && ~isD(f)), continue, end
    V3=fu(f,:);                           % the three global face vertices
    dofs=[];sgn=[1,-1];ok=true;
    for side=1:cnt(f)
        iel=f2el(f,side);j=f2loc(f,side);iv=tet(iel,:);
        oth=setdiff(1:4,j);
        loc=zeros(1,3);
        for m=1:3
            q=oth(iv(oth)==V3(m));
            if(isempty(q)), ok=false;break, end
            loc(m)=q;
        end
        if(~ok), break, end
        dofs=[dofs;fac(iel,loc)*sgn(side)];   % signed, one row per side
    end
    if(~ok), continue, end
    % d_X = sum over sides of sign * u_{face dof}
    ns=size(dofs,1);
    fd=zeros(1,3*ns);sg=zeros(1,3*ns);
    for side=1:ns
        fd((side-1)*3+(1:3))=abs(dofs(side,:));
        sg((side-1)*3+(1:3))=sign(dofs(side,1));
    end
    % assemble  w * sum_{X,Y} MF0(X,Y) d_X d_Y   for each component
    w=2*my*opts.gamma1*farea(f)/hF(f);
    n=3*ns;
    Mbig=zeros(n);
    for a=1:n
        for b=1:n
            Xa=mod(a-1,3)+1;Xb=mod(b-1,3)+1;
            Mbig(a,b)=w*sg(a)*sg(b)*MF0(Xa,Xb);
        end
    end
    for comp=0:2
        eqs=3*fd-2+comp;
        [I,J]=ndgrid(eqs,eqs);
        rp(np+1:np+n*n)=I(:);cp(np+1:np+n*n)=J(:);vp(np+1:np+n*n)=Mbig(:);
        np=np+n*n;
    end
end
K=sparse([ra(1:ua);rp(1:np)],[ca(1:ua);cp(1:np)],[va(1:ua);vp(1:np)],neqU,neqU);
B=sparse(rb(1:ub),cb(1:ub),vb(1:ub),nele,neqU);
M=spdiags(vol,0,nele,nele);Mi=spdiags(1./vol,0,nele,nele);

% ---------------- boundary conditions ----------------
if(nargin < 6 || isempty(dind))
    error('solvedisccrs3: prescribe dind and dval');
end
[ind,ia]=unique(dind(:),'first');bcval=dval(ia);
int=setdiff((1:neq)',ind);
ee=(1:nele)';

u=zeros(neqU,1);p=zeros(nele,1);act=true(nele,1);actold=act;
conv=false;nit=0;
for it=1:100
    if(it==1), act=true(nele,1); else, act=(gam*p-(B*u)./vol>=0); end
    if(it>1 && isequal(act,actold)), conv=true;break, end
    actold=act;
    EA=sparse(ee(act),ee(act),1,nele,nele);EI=sparse(ee(~act),ee(~act),1,nele,nele);
    Smat=[K+(1/gam)*B'*EA*Mi*B, -B'*EA; -EA*B, -gam*EI*M];
    fr=-Smat(:,ind)*bcval;
    v=zeros(neq,1);v(ind)=bcval;
    v(int)=Smat(int,int)\fr(int);
    u=v(1:neqU);p=v(neqU+1:end);nit=it;
end

ux=u(1:3:end);uy=u(2:3:end);uz=u(3:3:end);
D=(B*u)./vol;
checktol=1e-9*max([1,norm(p,inf),norm(D,inf)]);
signok=all(p>=-checktol) && all(D>=-checktol);
info=struct('iterations',nit,'converged',conv,'gamma',gam,'gamma1',opts.gamma1, ...
   'ncav',nnz(~act),'nele',nele,'cavfrac',nnz(~act)/nele,'act',act,'div',D, ...
   'faces',fu,'xf',xf,'yf',yf,'zf',zf,'onbnd',onbnd,'vol',vol,'farea',farea, ...
   'compl',max(abs(p.*D)),'signok',signok,'checktol',checktol, ...
   'pmin',min(p),'pmax',max(p),'ndof',neq,'nface',nface);
if(opts.verbose)
  fprintf(['solvedisccrs3: %d its (conv %d), gamma = %g, gamma1 = %g, %d dof\n' ...
     '   p in [%.4g, %.4g], min div u = %.3e, max|p*div u| = %.3e\n'], ...
     nit,conv,gam,opts.gamma1,neq,min(p),max(p),min(D),max(abs(p.*D)));
end
end
