% pitcompareends.pdf (Figure 6): pressure on the mid plane for the four pits
% of Table 8, Reynolds against Stokes, with pressure-normal-flow Stokes ends
% and the steepness increasing downwards.  The bars below each profile mark
% the cavitated samples; the shaded strips lie outside the window 4 < x < 20.
% Requires Optimization Toolbox (quadprog).
clear all, close all
here=fileparts(mfilename('fullpath'));addpath(here);
configure_submission_figures

cases=[0.5 8;1 4;1 2;1 1];L=24;c=1;
res=cell(size(cases,1),1);
for k=1:size(cases,1)
    res{k}=pressure_end_pit(cases(k,1),cases(k,2));
    fprintf('delta/r %.4f: peaks %.5f/%.5f, fronts %.4f/%.4f\n',cases(k,1)/cases(k,2), ...
        res{k}.peak_reynolds,res{k}.peak_stokes,res{k}.front_reynolds,res{k}.front_stokes);
end

figure(1),clf
set(gcf,'Units','centimeters','Position',[2 2 17 19],'Color','w')
for k=1:4
    R=res{k};
    ax=subplot(4,1,k);hold on,box on
    ym=max(R.pRf)*1.15;yl=[-0.15*max(R.pRf) ym];
    patch([0 4*c 4*c 0],[yl(1) yl(1) yl(2) yl(2)],[0.9 0.9 0.9], ...
        'EdgeColor','none','HandleVisibility','off')
    patch([L-4*c L L L-4*c],[yl(1) yl(1) yl(2) yl(2)],[0.9 0.9 0.9], ...
        'EdgeColor','none','HandleVisibility','off')
    plot(R.xf,zeros(size(R.xf)),'k:','HandleVisibility','off')
    hR=plot(R.xf,R.pRf,'k-','LineWidth',1.3);
    hS=plot(R.xf,R.pSf,'r--','LineWidth',1.1);
    % cavitated samples as line segments broken by NaN; markers would pull a
    % MathWorks glyph font into the PDF
    w=(R.xf>4*c)&(R.xf<L-4*c);yb1=yl(1)*0.45;yb2=yl(1)*0.80;
    br=nan(size(R.xf));br(w&R.cRf)=yb1;
    bs=nan(size(R.xf));bs(w&R.cSf)=yb2;
    plot(R.xf,br,'k-','LineWidth',2.2,'HandleVisibility','off')
    plot(R.xf,bs,'r-','LineWidth',2.2,'HandleVisibility','off')
    xlim([0 L]),ylim(yl),ylabel('{\itp}')
    if k==1
        legend([hR hS],{'Reynolds','Stokes'},'Location','NorthWest', ...
            'Box','off','FontSize',9);
    end
    if k<4, set(ax,'XTickLabel',[]); end
    set(ax,'FontSize',9,'Layer','top')
end
xlabel('{\itx}')
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

outdir=figure_output_dir(here);
exportgraphics(gcf,fullfile(outdir,'pitcompareends.pdf'), ...
    'ContentType','vector','BackgroundColor','white')
fprintf('wrote %s\n',fullfile(outdir,'pitcompareends.pdf'));
