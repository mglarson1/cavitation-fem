function results=reproduce_domain_length(outdir)
% Table 12: fixed-mesh, fixed-window inlet and outlet location experiment.
if nargin<1,outdir=fullfile(pwd,'results');end
folder=fullfile(outdir,'domain_length_2026-09-08');
verify_domain_length(folder,[0,24;-12,36;-36,60]);
verify_domain_length(fullfile(folder,'one_sided'),[0,24;-12,24;0,36]);
results=audit_domain_length_results(folder);
end
