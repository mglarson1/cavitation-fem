function results=audit_pressure_end_results(folder)
% Read-only numerical acceptance checks; write manuscript table summaries.
if nargin<1,folder=fullfile(pwd,'results','pressure_ends_2026-09-08');end
tags={'case_0_24_m1_g0_tol1e-11','case_-12_36_m1_g0_tol1e-11', ...
    'case_-36_60_m1_g0_tol1e-11','case_0_24_m2_g0_tol1e-11', ...
    'case_0_24_m4_g0_tol1e-11','case_0_24_m4_g1_tol1e-11', ...
    'case_0_24_m4_g0_tol1e-11_shift','case_0_24_m1_g0_tol1e-12'};
rows=cell(numel(tags),1);central=cell(3,1);grids=cell(numel(tags),1);
for j=1:numel(tags)
    s=load(fullfile(folder,[tags{j},'.mat']));
    assert(s.iR.converged && s.iR.signok && s.iS.converged && s.iS.signok);
    assert(s.iS.qp.qp_flag>0 && s.iS.qp.converged);
    tol=s.iS.qp.kkt_tolerance;
    assert(min(s.pS)>=-tol/100 && min(s.iS.div)>=-tol);
    assert(s.iS.stationarity<1e-10 && s.iS.compl<1e-10);
    assert(isequaln(s.row.stokes_front,s.row.stokes_lower_trace_front));
    assert(all(isfinite(s.profiles.stokes_pressure)));
    assert(s.row.reynolds_segments==1 && s.row.reynolds_left_censored);
    if s.row.stokes_length>0
        assert(s.row.stokes_segments==1 && s.row.stokes_left_censored);
        assert(abs(s.row.stokes_front_upper-s.row.stokes_front-.000375)<1e-12);
    end
    if ~ismember('shift',s.row.Properties.VariableNames),s.row.shift=false;end
    rows{j}=s.row;grids{j}=s.xx;
    % Domain extension must preserve every central triangle coordinate.
    if j<=3
        xx=s.x(s.tri);zz=s.z(s.tri);
        keep=all(xx>=4 & xx<=20,2);
        central{j}=sortrows([xx(keep,:),zz(keep,:)]);
    end
end
assert(isequal(central{1},central{2}) && isequal(central{1},central{3}));
assert(isequal(grids{5}(grids{5}>=6 & grids{5}<=18), ...
               grids{6}(grids{6}>=6 & grids{6}<=18)));
assert(abs(rows{1}.stokes_front-rows{8}.stokes_front)<1e-12);
assert(abs(rows{5}.stokes_front-rows{6}.stokes_front)<1e-12);
assert(abs(rows{5}.stokes_peak-rows{6}.stokes_peak)<3e-5);
assert(abs(rows{7}.stokes_front-rows{5}.stokes_front)>.01);
assert(rows{3}.stokes_length==0);
s=load(fullfile(folder,[tags{3},'.mat']),'profiles');
assert(min(s.profiles.stokes_pressure)>4e-3);
flat=readtable(fullfile(folder,'flat_channel.csv'));
assert(height(flat)==4);
assert(max(flat{:,{'ux_error','uz_error','pressure_error'}},[],'all')<1e-12);
assert(max(flat.exact_stationarity)<4e-15 && max(flat.exact_divergence)==0);
refs=readtable(fullfile(folder,'reynolds_front_reference.csv'));
assert(max(abs(refs.closure_residual))<1e-12);
ref=refs.front_relative_to_pit(refs.right==24);
assert(abs(rows{5}.reynolds_front-ref)<.01);
assert(abs(rows{7}.reynolds_front-ref)<.0025);
results=struct('domain',vertcat(rows{1:3}), ...
    'refinement',vertcat(rows{[1,4,5,6,7]}),'all',vertcat(rows{:}));
writetable(results.domain,fullfile(folder,'table13_pressure_ends.csv'));
writetable(results.refinement,fullfile(folder,'table14_front_refinement.csv'));
writetable(results.all,fullfile(folder,'accepted_cases.csv'));
fprintf('Pressure-end audit passed: four flat checks and eight pit solves.\n');
end
