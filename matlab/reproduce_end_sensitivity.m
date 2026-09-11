function table10 = reproduce_end_sensitivity(outdir)
% Reproduce manuscript Table 10: end sensitivity of the steepest pit on five
% domains with identical meshes on [0,24].  The traction-free Stokes columns
% are read from the audit written by reproduce_domain_length, which must be
% run first; the pressure-normal-flow columns are computed here.  Requires
% Optimization Toolbox (quadprog).

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
tfile=fullfile(outdir,'domain_length_2026-09-08','table12_domain_length.csv');
if(~isfile(tfile))
    error('reproduce:missing','Run reproduce_domain_length first (%s).',tfile)
end
tf=readtable(tfile);
folder=fullfile(outdir,'end_sensitivity');
D=[0 24;-12 24;0 36;-12 36;-36 60];n=size(D,1);
left=D(:,1);right=D(:,2);
[reynolds_peak,reynolds_front,traction_free_stokes_peak,traction_free_stokes_cavity, ...
    pressure_normal_flow_stokes_peak,pressure_normal_flow_stokes_front]=deal(nan(n,1));
for j=1:n
    row=run_pressure_end_case(struct('left',D(j,1),'right',D(j,2)),folder);
    reynolds_peak(j)=row.reynolds_peak;reynolds_front(j)=row.reynolds_front;
    pressure_normal_flow_stokes_peak(j)=row.stokes_peak;
    pressure_normal_flow_stokes_front(j)=row.stokes_front;
    i=find(tf.left==D(j,1)&tf.right==D(j,2));assert(numel(i)==1)
    traction_free_stokes_peak(j)=tf.stokes_peak(i);
    traction_free_stokes_cavity(j)=tf.stokes_cavity(i);
end
prof=readtable(fullfile(folder,'case_-36_60_m1_g0_tol1e-11_profiles.csv'));
min_window_pressure_longest=min(prof.stokes_pressure);
table10=table(left,right,reynolds_peak,reynolds_front,traction_free_stokes_peak, ...
    traction_free_stokes_cavity,pressure_normal_flow_stokes_peak,pressure_normal_flow_stokes_front);
writetable(table10,fullfile(outdir,'end_sensitivity.csv'));

check_close('Table 10 Reynolds peaks',reynolds_peak,[2.2158;2.2158;2.4018;2.4018;2.5115],6e-5,0);
check_close('Table 10 Reynolds fronts',reynolds_front,[-1.781;-1.781;-1.969;-1.969;-2.156],6e-4,0);
check_close('Table 10 traction-free Stokes peaks',traction_free_stokes_peak, ...
    [2.1473;2.1472;2.3488;2.3474;2.4659],6e-5,0);
check_close('Table 10 traction-free Stokes cavities',traction_free_stokes_cavity, ...
    [0.843;1.218;0;0;0],2e-4,0);
check_close('Table 10 pressure-normal-flow Stokes peaks',pressure_normal_flow_stokes_peak, ...
    [2.1289;2.1289;2.3353;2.3353;2.4561],6e-5,0);
check_close('Table 10 pressure-normal-flow Stokes fronts',pressure_normal_flow_stokes_front(1:4), ...
    [-1.969;-1.969;-2.438;-2.438],6e-4,0);
assert(isnan(pressure_normal_flow_stokes_front(5)))
assert(abs(min_window_pressure_longest-4.6e-3)<1e-4)
fprintf('Reproduced Table 10 in %s\n',outdir)
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
