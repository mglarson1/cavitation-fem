function results=reproduce_pressure_ends(outdir)
% Tables 13--14: Couette calibration, compatible ends and front refinement.
% Requires Optimization Toolbox. The full run takes about 8 minutes.
if nargin<1,outdir=fullfile(pwd,'results');end
folder=fullfile(outdir,'pressure_ends_2026-09-08');
verify_pressure_end_flat_channel(folder);
refs=cell(3,1);domains=[0,24;-12,36;-36,60];
for j=1:3
    run_pressure_end_case(struct('left',domains(j,1),'right',domains(j,2)),folder);
    refs{j}=reynolds_pit_front_reference(domains(j,2));
end
writetable(vertcat(refs{:}),fullfile(folder,'reynolds_front_reference.csv'));
for m=[2,4]
    run_pressure_end_case(struct('left',0,'right',24,'factor',m),folder);
end
run_pressure_end_case(struct('left',0,'right',24,'factor',4,'global_refine',true),folder);
run_pressure_end_case(struct('left',0,'right',24,'factor',4,'shift',true),folder);
run_pressure_end_case(struct('left',0,'right',24,'kkt_tolerance',1e-12),folder);
results=audit_pressure_end_results(folder);
end
