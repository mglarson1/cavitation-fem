function reproduce(mode)
% Reproduce numerical tables and figures from one public entry point.
%
%   reproduce('quick')    Tables 1--4
%   reproduce('tables')   Tables 1--12, the 48 smoothed Newton runs and the
%                         mechanical-pressure values of Section 5.6
%   reproduce('figures')  All seven manuscript figures
%   reproduce('korn')     Korn constants of Remark A.6 (about 6 GB memory)
%   reproduce('full')     Tables, figures and Korn constants
%   reproduce('verification') 2D Reynolds, unsmoothed QP and end studies
%   reproduce('domain')   Traction-free and calibrated end sensitivity (Table 10)
%   reproduce('pressure-ends') Couette calibration and front refinement (Table 9)
%   reproduce('certify')  Tables, Korn constants and all prose checks
%
% The default is 'quick'.  All CSV files are written to matlab/results.
% Figures go to tex/figures in the manuscript tree and to figures in a
% standalone code release.  Set CAVITATION_FIGURE_DIR to override this.
% Optimization Toolbox (quadprog) is required for Tables 8--10 and 12.

if(nargin<1 || isempty(mode)), mode='quick'; end
mode=lower(string(mode));
here=fileparts(mfilename('fullpath'));addpath(here);
outdir=fullfile(here,'results');

if(any(mode==["quick","tables","full","certify"]))
    fprintf('\n=== Tables 1--3: Reynolds verification ===\n')
    reproduce_reynolds_tables(outdir);
    fprintf('\n=== Table 4: two-dimensional Stokes channel ===\n')
    reproduce_stokes_core_tables(outdir);
end

if(any(mode==["tables","full","certify"]))
    fprintf('\n=== Table 5: three-dimensional verification ===\n')
    reproduce_stokes_3d_table(outdir);
    fprintf('\n=== Table 6 and mechanical pressure: constitutive laws ===\n')
    reproduce_constitutive_table(outdir);
    verify_mechanical_pressure(outdir);
    fprintf('\n=== Table 7: discretizations ===\n')
    reproduce_discretization_table(outdir);
    fprintf('\n=== Table 8: Reynolds--Stokes comparison ===\n')
    reproduce_model_comparison(outdir);
    fprintf('\n=== Table 11: active set with sign rule and smoothed Newton ===\n')
    reproduce_solver_table(outdir);
    fprintf('\n=== 48 smoothed Newton runs ===\n')
    reproduce_smoothed_newton_table(outdir);
end

if(any(mode==["tables","full","certify","verification"]))
    if(mode=="verification"), verify_reynolds_2d_convergence(outdir); end
    fprintf('\n=== Unsmoothed Stokes QP reference ===\n')
    verify_stokes_qp_reference(outdir);
    fprintf('\n=== Table 12: smoothing against the QP reference ===\n')
    verify_smoothing_against_qp(outdir);
    verify_smoothing_against_qp(outdir,[1e-16,1e-18],[1e-10,1e-12], ...
        'smoothing_cavity_tolerance.csv');
end

if(any(mode==["pressure-ends","tables","full","certify","verification"]))
    fprintf('\n=== Table 9: Couette calibration and cavity fronts ===\n')
    reproduce_pressure_ends(outdir);
end

if(any(mode==["domain","tables","full","certify","verification"]))
    fprintf('\n=== Table 10: end sensitivity ===\n')
    reproduce_domain_length(outdir);
    reproduce_end_sensitivity(outdir);
end

if(any(mode==["korn","full","certify"]))
    fprintf('\n=== Remark A.6: Korn constants ===\n')
    verify_korn_constants(outdir);
end

if(mode=="certify")
    fprintf('\n=== Quantitative claims stated in prose ===\n')
    verify_manuscript_claims(outdir,'full');
    fprintf('\n=== Safeguarded Newton matrix ===\n')
    verify_newton_safeguard(outdir);
end

if(any(mode==["figures","full"]))
    fprintf('\n=== Figures 1--3: Reynolds pit ===\n')
    makepit
    makepitlambda
    makepitline
    fprintf('\n=== Figure 4: Stokes channel ===\n')
    makechannel
    fprintf('\n=== Figure 5: three-dimensional Stokes verification ===\n')
    makestokes3d
    fprintf('\n=== Figures 6--7: Reynolds--Stokes pit comparison ===\n')
    makepitcompare
    makepitfield
end
