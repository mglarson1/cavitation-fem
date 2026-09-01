function results = reproduce_stokes_core_tables(outdir)
% Reproduce manuscript Tables 5--6 and the steep-pit mesh sensitivity.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

%% Table 5: two-dimensional channel refinement
refs=(2:5)';n=numel(refs);
elements=zeros(n,1);ndof=elements;iterations=elements;pmax=elements;
cavfrac=elements;complementarity=elements;
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');
for k=1:n
    nx=3*2^refs(k);ny=2^refs(k);
    [tri,x,y]=meshrect(3,1,nx,ny);
    [dind,dval]=channel_cr_boundary(tri,x,y);
    [~,~,p,info]=solvedisccrs(tri,x,y,100,dind,dval,opCR);

    % The P0 pressure has two traces on y=1/2.  The paper uses the upper
    % trace at longitudinal cell centers, avoiding points on mesh edges.
    hx=3/nx;xs=((0:nx-1)'+0.5)*hx;
    TR=triangulation(tri,x,y);
    eid=pointLocation(TR,[xs,(0.5+1e-10)*ones(size(xs))]);
    pmax(k)=max(p(eid(xs<2.5)));
    elements(k)=size(tri,1);ndof(k)=info.ndof;
    iterations(k)=info.iterations;cavfrac(k)=info.cavfrac;
    complementarity(k)=info.compl;
    assert(info.converged && info.signok)
end
table05=table(refs,elements,ndof,iterations,pmax,cavfrac,complementarity);
writetable(table05,fullfile(outdir,'table05_stokes_channel_refinement.csv'));
check_close('Table 5 elements',elements,[96;384;1536;6144],0,0);
check_close('Table 5 degrees of freedom',ndof,[416;1600;6272;24832],0,0);
check_close('Table 5 iterations',iterations,[3;2;4;3],0,0);
check_close('Table 5 peak pressures',pmax, ...
    [5.72708677;5.77392315;5.82275035;5.85246269],6e-8,0);
check_close('Table 5 cavity fractions',cavfrac, ...
    [0.0104166667;0.0208333333;0.0188802083;0.0187174479],6e-10,0);
assert(max(complementarity)<5e-13)

%% Table 6: Reynolds--Stokes comparison on four common pits
cases=[0.5 8;1 4;1 2;1 1];n=size(cases,1);
delta=cases(:,1);radius=cases(:,2);steepness=delta./radius;
peak_reynolds=zeros(n,1);peak_stokes=peak_reynolds;profile_difference_pct=peak_reynolds;
cavity_reynolds=peak_reynolds;cavity_stokes=peak_reynolds;
for k=1:n
    row=common_pit_metrics(delta(k),radius(k),256,16,'cr',opCR);
    peak_reynolds(k)=row.peakR;peak_stokes(k)=row.peakS;
    profile_difference_pct(k)=row.profileDiff;
    cavity_reynolds(k)=row.cavLenR;cavity_stokes(k)=row.cavLenS;
end
table06=table(delta,radius,steepness,peak_reynolds,peak_stokes, ...
    profile_difference_pct,cavity_reynolds,cavity_stokes);
writetable(table06,fullfile(outdir,'table06_reynolds_stokes_comparison.csv'));
check_close('Table 6 Reynolds peaks',peak_reynolds, ...
    [3.3653;5.4940;3.8098;2.2158],6e-5,0);
check_close('Table 6 Stokes peaks',peak_stokes, ...
    [3.3730;5.4583;3.7572;2.1473],6e-5,0);
check_close('Table 6 profile differences',profile_difference_pct, ...
    [0.80;1.30;2.26;4.10],6e-3,0);
check_close('Table 6 Reynolds cavity lengths',cavity_reynolds, ...
    [1.905;2.748;4.812;6.219],2e-4,0);
check_close('Table 6 Stokes cavity lengths',cavity_stokes, ...
    [2.187;2.937;2.061;0.843],2e-4,0);

%% Steep-pit sensitivity quoted immediately after Table 6
nxvals=[128;256;512];nzvals=[8;16;32];n=numel(nxvals);
elements=zeros(n,1);peak_stokes=zeros(n,1);cavity_stokes=zeros(n,1);
ratio=zeros(n,1);
for k=1:n
    row=common_pit_metrics(1,1,nxvals(k),nzvals(k),'cr',opCR);
    elements(k)=2*nxvals(k)*nzvals(k);peak_stokes(k)=row.peakS;
    cavity_stokes(k)=row.cavLenS;
end
% As stated in the manuscript, compare every Stokes mesh with the
% Table-6 Reynolds length on the fixed 256-interval grid.
ratio=cavity_reynolds(end)./cavity_stokes;
table06_sensitivity=table(nxvals,nzvals,elements,peak_stokes,cavity_stokes,ratio);
writetable(table06_sensitivity, ...
    fullfile(outdir,'table06_steep_pit_mesh_sensitivity.csv'));
check_close('Table 6 sensitivity cavity lengths',cavity_stokes, ...
    [0.750;0.843;0.798],2e-4,0);
check_close('Table 6 sensitivity ratios',ratio,[8.292;7.377;7.793],6e-4,0);

results=struct('table05',table05,'table06',table06, ...
    'table06_sensitivity',table06_sensitivity);
fprintf('Reproduced Tables 5--6 in %s\n',outdir)
end

function [dind,dval]=channel_cr_boundary(tri,x,y)
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
ya=y(eu(:,1));yb=y(eu(:,2));xa=x(eu(:,1));xb=x(eu(:,2));
gmean=(ya+yb)/2-(ya.^2+ya.*yb+yb.^2)/3;
wall=find(onb&((abs(ya)<1e-12&abs(yb)<1e-12)| ...
    (abs(ya-1)<1e-12&abs(yb-1)<1e-12)));
inflow=find(onb&abs(xa)<1e-12&abs(xb)<1e-12);
dind=[2*wall-1;2*wall;2*inflow-1;2*inflow];
dval=[zeros(2*numel(wall),1);gmean(inflow);zeros(numel(inflow),1)];
end

function row=common_pit_metrics(delta,rp,nx,nz,method,opts)
L=24;xp=L/2;c=1;V=1;gam=100;ny=4;
Hf=@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
dfun=@(x,y) Hf(x)/c;
ffun=@(x,y) -(-2*(x-xp)/rp^2).*delta.*exp(-((x-xp)/rp).^2);

% One-dimensional Reynolds solution represented on a narrow strip.
[triR,xR,yR]=meshrect(L,ny*L/nx,nx,ny);
dnR=find(abs(xR)<1e-9|abs(xR-L)<1e-9);
evalc('[P,~,iR]=solvereynolds(triR,xR,yR,dfun,ffun,1,dnR);');
pR=6*P;mid=find(abs(yR-(ny*L/nx)/2)<1e-9);
[xRs,is]=sort(xR(mid));mid=mid(is);pRs=pR(mid);actRs=iR.act(mid);

% Stokes cross section.  The flat upper surface slides and the shaped
% lower surface is stationary; both ends are traction free.
[triS,xh,zh]=meshrect(L,1,nx,nz);xm=xh;zm=c-zh.*Hf(xh);
E=[triS(:,[2,3]);triS(:,[1,3]);triS(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
bot=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
top=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*bot-1;2*bot;2*top-1;2*top];
dval=[V*ones(numel(bot),1);zeros(numel(bot),1); ...
    zeros(numel(top),1);zeros(numel(top),1)];
switch method
    case 'cr'
        [~,~,pS,iS]=solvedisccrs(triS,xm,zm,gam,dind,dval,opts);
    case 'cr-unstabilized'
        evalc('[~,~,pS,iS]=solvedisccr(triS,xm,zm,gam,dind,dval);');
    otherwise
        error('Unknown pit method %s',method)
end

xs=linspace(0,L,8001)';TRm=triangulation(triS,xm,zm);
eid=pointLocation(TRm,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
pSs=nan(size(xs));pSs(ok)=pS(eid(ok));
cavS=false(size(xs));cavS(ok)=~iS.act(eid(ok));
pRi=interp1(xRs,pRs,xs);
ibin=discretize(xs,xRs);valid=~isnan(ibin);cavR=false(size(xs));
cavR(valid)=actRs(ibin(valid))&actRs(ibin(valid)+1);
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);
row=struct('peakR',max(pRi(w)),'peakS',max(pSs(w)), ...
    'profileDiff',100*max(abs(pSs(w)-pRi(w)))/max(pRi(w)), ...
    'cavLenR',dx*nnz(cavR(w)),'cavLenS',dx*nnz(cavS(w)), ...
    'infoS',iS);
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
