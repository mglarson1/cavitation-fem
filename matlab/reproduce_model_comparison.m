function table08 = reproduce_model_comparison(outdir)
% Reproduce manuscript Table 8: Reynolds--Stokes comparison for the four
% pits with pressure-normal-flow Stokes ends on [0,24].  Requires
% Optimization Toolbox (quadprog).

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

cases=[0.5 8;1 4;1 2;1 1];n=size(cases,1);
delta=cases(:,1);radius=cases(:,2);steepness=delta./radius;
[peak_reynolds,peak_stokes,profile_difference_pct,front_reynolds, ...
    front_stokes,min_pressure,min_divergence,complementarity]=deal(zeros(n,1));
left_censored=false(n,1);
for k=1:n
    c=pressure_end_pit(delta(k),radius(k));
    peak_reynolds(k)=c.peak_reynolds;peak_stokes(k)=c.peak_stokes;
    profile_difference_pct(k)=c.profile_difference_pct;
    front_reynolds(k)=c.front_reynolds;front_stokes(k)=c.front_stokes;
    left_censored(k)=c.left_censored;min_pressure(k)=c.min_pressure;
    min_divergence(k)=c.min_divergence;complementarity(k)=c.complementarity;
end
table08=table(delta,radius,steepness,peak_reynolds,peak_stokes, ...
    profile_difference_pct,front_reynolds,front_stokes,left_censored, ...
    min_pressure,min_divergence,complementarity);
writetable(table08,fullfile(outdir,'model_comparison_pressure_normal_flow.csv'));

check_close('Table 8 Reynolds peaks',peak_reynolds,[3.3653;5.4940;3.8098;2.2158],6e-5,0);
check_close('Table 8 Stokes peaks',peak_stokes,[3.3631;5.4465;3.7366;2.1289],6e-5,0);
check_close('Table 8 profile differences',profile_difference_pct,[0.63;1.34;2.48;4.52],6e-3,0);
check_close('Table 8 Reynolds fronts',front_reynolds,[-6.094;-5.250;-3.188;-1.781],6e-4,0);
check_close('Table 8 Stokes fronts',front_stokes,[-5.813;-5.063;-3.094;-1.969],6e-4,0);
assert(all(left_censored))
assert(min(min_pressure)>=-1e-11 && min(min_divergence)>=-1e-11 && max(complementarity)<1e-12)
fprintf('Reproduced Table 8 in %s\n',outdir)
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
