function [ux,uy,uz,p,info]=solvedisccrs3_ip(tet,xnod,ynod,znod,gam,dind,dval,opts)
% solvedisccrs3.m with a memory-lean inner solver, for fine meshes (used for
% the 36x12x12 mesh of Table 5 and Figure 5).  The
% discretization, the active set iteration and the output are those of
% solvedisccrs3.m.  The only change is how the frozen system
%       K_A u - B_A' p = f ,   B_A u = g ,   p = 0 on cavitated elements,
% is solved.  Instead of a sparse LDL factorization of the indefinite saddle
% matrix, the iterated penalty (augmented Lagrangian) method
%       (K_A + r B_A' M^-1 B_A) u^{k+1} = f + B_A' p^k + r B_A' M^-1 g ,
%       p^{k+1} = p^k - r M^-1 (B_A u^{k+1} - g) ,
% is used, with one sparse Cholesky factorization per active set.  The
% fixed point is the solution of the frozen system, and since div V^h = Q^h
% the pressure error contracts by a factor O(1/r) per step.
%
%   extra opts: .r      penalty (default 1e4)
%               .ptol   stop when max|p^{k+1}-p^k| <= ptol*max(1,|p|_inf)
%                       (default 1e-10), or when the update stagnates
%               .order  'dissect' (nested dissection, default) or 'amd'

if(nargin < 5 || isempty(gam)), gam=100; end
if(nargin < 8), opts=struct; end
if(~isfield(opts,'gamma1')),  opts.gamma1=1;       end
if(~isfield(opts,'lamfac')),  opts.lamfac=-2/3;    end
if(~isfield(opts,'verbose')), opts.verbose=1;      end
if(~isfield(opts,'r')),       opts.r=1e4;          end
if(~isfield(opts,'ptol')),    opts.ptol=1e-10;     end
if(~isfield(opts,'order')),   opts.order='dissect';end
tstart=tic;

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
clear F

xf=(xnod(fu(:,1))+xnod(fu(:,2))+xnod(fu(:,3)))/3;
yf=(ynod(fu(:,1))+ynod(fu(:,2))+ynod(fu(:,3)))/3;
zf=(znod(fu(:,1))+znod(fu(:,2))+znod(fu(:,3)))/3;

v1=[xnod(fu(:,2))-xnod(fu(:,1)),ynod(fu(:,2))-ynod(fu(:,1)),znod(fu(:,2))-znod(fu(:,1))];
v2=[xnod(fu(:,3))-xnod(fu(:,1)),ynod(fu(:,3))-ynod(fu(:,1)),znod(fu(:,3))-znod(fu(:,1))];
cr=cross(v1,v2,2);
farea=0.5*sqrt(sum(cr.^2,2));
clear v1 v2 cr

% ---------------- elementwise forms ----------------
[ra,ca,va,ua]=assemble(nele*144);
rb=zeros(nele*12,1);cb=rb;vb=rb;ub=0;vol=zeros(nele,1);
for iel=1:nele
    iv=tet(iel,:);xc=xnod(iv);yc=ynod(iv);zc=znod(iv);
    [~,lx,ly,lz,V]=basis3(sum(xc)/4,sum(yc)/4,sum(zc)/4,xc,yc,zc);
    px=-3*lx;py=-3*ly;pz=-3*lz;
    vol(iel)=V;
    Be=zeros(6,12);
    Be(1,1:3:end)=px';
    Be(2,2:3:end)=py';
    Be(3,3:3:end)=pz';
    Be(4,1:3:end)=py';Be(4,2:3:end)=px';
    Be(5,1:3:end)=pz';Be(5,3:3:end)=px';
    Be(6,2:3:end)=pz';Be(6,3:3:end)=py';
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

MF0=(3*eye(3)-ones(3))/4;
np=0;maxn=nface*3*36;                     % at most 6 x 6 entries per component
rp=zeros(maxn,1);cp=rp;vp=rp;
for f=1:nface
    if(onbnd(f) && ~isD(f)), continue, end
    V3=fu(f,:);
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
        dofs=[dofs;fac(iel,loc)*sgn(side)];
    end
    if(~ok), continue, end
    ns=size(dofs,1);
    fd=zeros(1,3*ns);sg=zeros(1,3*ns);
    for side=1:ns
        fd((side-1)*3+(1:3))=abs(dofs(side,:));
        sg((side-1)*3+(1:3))=sign(dofs(side,1));
    end
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
clear ra ca va rp cp vp rb cb vb f2el f2loc
tass=toc(tstart);

% ---------------- boundary conditions ----------------
if(nargin < 6 || isempty(dind))
    error('solvedisccrs3_ip: prescribe dind and dval');
end
[ind,ia]=unique(dind(:),'first');bcval=dval(ia);
fre=setdiff((1:neqU)',ind);nfre=numel(fre);
Kff=K(fre,fre);Kfd=K(fre,ind);
Bf=B(:,fre);gB=-B(:,ind)*bcval;          % B*u = 0  <=>  Bf*u_f = gB
r=opts.r;

u=zeros(neqU,1);u(ind)=bcval;p=zeros(nele,1);act=true(nele,1);actold=act;
conv=false;nit=0;inner=zeros(0,1);nnzL=0;tfac=0;lastres=NaN;lastdp=NaN;
for it=1:100
    if(it==1), act=true(nele,1); else, act=(gam*p-(B*u)./vol>=0); end
    if(it>1 && isequal(act,actold)), conv=true;break, end
    actold=act;
    BA=Bf(act,:);gA=gB(act);wA=1./vol(act);
    BtW=BA'*spdiags(wA,0,numel(wA),numel(wA));
    BWB=BtW*BA;
    fA=-Kfd*bcval+(1/gam)*(BtW*gA);
    Ar=Kff+(1/gam+r)*BWB;
    clear BWB
    L=[];t0=tic;
    switch opts.order
      case 'dissect', q=dissect(Ar);[L,flag]=chol(Ar(q,q),'lower');
      otherwise,      [L,flag,q]=chol(Ar,'lower','vector');
    end
    if(flag), error('solvedisccrs3_ip: Cholesky failed (flag %d)',flag); end
    tfac=tfac+toc(t0);nnzL=max(nnzL,nnz(L));
    clear Ar
    pA=p(act);                             % warm start
    c0=fA+r*(BtW*gA);
    dpold=inf;
    for k=1:50
        rhs=c0+BA'*pA;
        uf=zeros(nfre,1);uf(q)=L'\(L\rhs(q));
        res=(BA*uf-gA).*wA;                % divergence on active elements
        dp=-r*res;pA=pA+dp;
        ndp=norm(dp,inf);
        if(ndp<=opts.ptol*max(1,norm(pA,inf)) || (k>2 && ndp>0.5*dpold)), break, end
        dpold=ndp;
    end
    u(fre)=uf;p=zeros(nele,1);p(act)=pA;nit=it;inner(end+1,1)=k;
    lastres=norm(res,inf);lastdp=ndp;
    if(opts.verbose)
      fprintf('  active set %d: %d cavitated, nnz(L) = %.4g, %d penalty its, max|div| = %.1e, max|dp| = %.1e\n', ...
          it,nnz(~act),nnz(L),k,lastres,lastdp);
    end
end
L=[];

ux=u(1:3:end);uy=u(2:3:end);uz=u(3:3:end);
D=(B*u)./vol;
checktol=1e-9*max([1,norm(p,inf),norm(D,inf)]);
signok=all(p>=-checktol) && all(D>=-checktol);
info=struct('iterations',nit,'converged',conv,'gamma',gam,'gamma1',opts.gamma1, ...
   'ncav',nnz(~act),'nele',nele,'cavfrac',nnz(~act)/nele,'act',act,'div',D, ...
   'faces',fu,'xf',xf,'yf',yf,'zf',zf,'onbnd',onbnd,'vol',vol,'farea',farea, ...
   'compl',max(abs(p.*D)),'signok',signok,'checktol',checktol, ...
   'pmin',min(p),'pmax',max(p),'ndof',neq,'nface',nface, ...
   'r',r,'inner',inner,'nnzL',nnzL,'tassemble',tass,'tfactor',tfac, ...
   'lastres',lastres,'lastdp',lastdp,'ttotal',toc(tstart));
if(opts.verbose)
  fprintf(['solvedisccrs3_ip: %d its (conv %d), gamma = %g, gamma1 = %g, r = %g, %d dof\n' ...
     '   p in [%.4g, %.4g], min div u = %.3e, max|p*div u| = %.3e\n' ...
     '   assembly %.1f s, factorizations %.1f s, total %.1f s\n'], ...
     nit,conv,gam,opts.gamma1,r,neq,min(p),max(p),min(D),max(abs(p.*D)), ...
     tass,tfac,info.ttotal);
end
end
