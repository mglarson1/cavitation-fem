% pitcompare.pdf: pressure on the mid plane for the four pits of Table
% tabC, Reynolds against Stokes, with the steepness increasing downwards.
% The Stokes computation is the stabilized Crouzeix-Raviart one of the
% paper; the figure block of mainpit.m, brought up to that discretization
% and with the excluded end regions shaded.
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

L=24;xp=L/2;c=1;V=1;gam=100;
nx=256;nz=16;          % Stokes cross section, as in Table tabC
ny=4;                  % transverse elements for the Reynolds strip
cases=[0.5 8; 1 4; 1 2; 1 1];
opCR=struct('gamma1',1,'lamfac',-2/3,'verbose',0,'solver','as');

res=cell(size(cases,1),1);
fprintf(['delta/r | max p: Reynolds  Stokes | profile diff. | ' ...
    'cavity length: Reynolds  Stokes | ratio\n'])
for k=1:size(cases,1)
    delta=cases(k,1);rp=cases(k,2);
    Hf =@(x) c*(1+delta*exp(-((x-xp)/rp).^2));
    dfun=@(x,y) Hf(x)/c;
    ffun=@(x,y) -(-2*(x-xp)/rp^2).*delta.*exp(-((x-xp)/rp).^2);

    % ---------- Reynolds, on a strip with P prescribed at the ends only
    [triR,xR,yR]=meshrect(L,ny*L/nx,nx,ny);
    dnR=find(abs(xR)<1e-9 | abs(xR-L)<1e-9);
    evalc('[P,lamR,iR]=solvereynolds(triR,xR,yR,dfun,ffun,1,dnR);');
    pR=6*P;
    mid=find(abs(yR-(ny*L/nx)/2)<1e-9);[xRs,is]=sort(xR(mid));
    mid=mid(is);pRs=pR(mid);actRs=iR.act(mid);

    % ---------- Stokes in the cross section, stabilized Crouzeix-Raviart
    % The pit is a depression in the stationary surface below; the flat
    % plane above slides.  zh = 0 is that plane, zh = 1 the shaped surface,
    % so the boundary conditions are assigned exactly as before and only the
    % geometric mapping is mirrored.  The problem is the reflection
    % z -> c-z of the other orientation and gives the same pressures.
    [triS,xh,zh]=meshrect(L,1,nx,nz);
    xm=xh; zm=c-zh.*Hf(xh);
    E=[triS(:,[2,3]);triS(:,[1,3]);triS(:,[1,2])];
    [eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
    bot=find(onb & abs(zh(eu(:,1)))<1e-9 & abs(zh(eu(:,2)))<1e-9);
    top=find(onb & abs(zh(eu(:,1))-1)<1e-9 & abs(zh(eu(:,2))-1)<1e-9);
    dind=[2*bot-1;2*bot;2*top-1;2*top];
    dval=[V*ones(numel(bot),1);zeros(numel(bot),1); ...
          zeros(numel(top),1);zeros(numel(top),1)];
    [~,~,pS,iS]=solvedisccrs(triS,xm,zm,gam,dind,dval,opCR);

    % sample the elementwise constant pressure on the mid plane
    TRm=triangulation(triS,xm,zm);
    xs=linspace(0,L,8001)';
    eid=pointLocation(TRm,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
    pSs=nan(size(xs));pSs(ok)=pS(eid(ok));
    cavS=false(size(xs));cavS(ok)=~iS.act(eid(ok));
    pRi=interp1(xRs,pRs,xs);
    ibin=discretize(xs,xRs);valid=~isnan(ibin);
    cavR=false(size(xs));
    cavR(valid)=actRs(ibin(valid)) & actRs(ibin(valid)+1);

    w=xs>4*c & xs<L-4*c;dx=xs(2)-xs(1);
    peakR=max(pRi(w));peakS=max(pSs(w));
    profileDiff=100*max(abs(pSs(w)-pRi(w)))/peakR;
    cavLenR=dx*nnz(cavR(w));cavLenS=dx*nnz(cavS(w));
    res{k}=struct('delta',delta,'r',rp,'xs',xs,'pS',pSs,'pR',pRi, ...
        'cavR',cavR,'cavS',cavS,'peakR',peakR,'peakS',peakS, ...
        'profileDiff',profileDiff,'cavLenR',cavLenR,'cavLenS',cavLenS);
    fprintf('%7.3f | %15.4f %7.4f | %12.2f%% | %23.3f %7.3f | %5.2f\n', ...
        delta/rp,peakR,peakS,profileDiff,cavLenR,cavLenS,cavLenR/cavLenS);
end

% ---------------- the figure ----------------
figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 17 19],'Color','w')
for k=1:4
    R=res{k};
    ax=subplot(4,1,k);hold on,box on
    ym=max(R.pR)*1.15;yl=[-0.15*max(R.pR) ym];
    % end regions excluded from the comparison
    patch([0 4*c 4*c 0],[yl(1) yl(1) yl(2) yl(2)],[0.9 0.9 0.9], ...
        'EdgeColor','none','HandleVisibility','off')
    patch([L-4*c L L L-4*c],[yl(1) yl(1) yl(2) yl(2)],[0.9 0.9 0.9], ...
        'EdgeColor','none','HandleVisibility','off')
    plot(R.xs,zeros(size(R.xs)),'k:','HandleVisibility','off')
    hR=plot(R.xs,R.pR,'k-','LineWidth',1.3);
    hS=plot(R.xs,R.pS,'r--','LineWidth',1.1);
    % the cavitated sets, which the pressure curves cannot show apart:
    % there the two models differ by "exactly zero" against "small but
    % positive", not by anything visible on this scale
    w=(R.xs>4*c)&(R.xs<L-4*c);
    yb1=yl(1)*0.45;yb2=yl(1)*0.80;
    % drawn as line segments, broken by NaN where the set is not cavitated,
    % rather than as dot markers: markers pull a MathWorks glyph font into
    % the PDF, and the segments are what the eye should read anyway
    br=nan(size(R.xs));br(w&R.cavR)=yb1;
    bs=nan(size(R.xs));bs(w&R.cavS)=yb2;
    plot(R.xs,br,'k-','LineWidth',2.2,'HandleVisibility','off')
    plot(R.xs,bs,'r-','LineWidth',2.2,'HandleVisibility','off')
    xlim([0 L]),ylim(yl)
    ylabel('{\itp}')
    if k==1
        legend([hR hS],{'Reynolds','Stokes'},'Location','NorthWest', ...
            'Box','off','FontSize',9);
    end
    if k<4, set(ax,'XTickLabel',[]); end
    set(ax,'FontSize',9,'Layer','top')
end
xlabel('{\itx}')
% all text in the figure in the body font of the paper, which elsarticle
% takes from txfonts, so that axis labels and captions do not clash
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pitcompare.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pitcompare.pdf'));
