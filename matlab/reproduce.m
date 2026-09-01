function reproduce(mode)
% Reproduce numerical tables and figures from one public entry point.
%
%   reproduce('quick')    Tables 1--6 and the quoted mesh sensitivity
%   reproduce('tables')   Tables 1--10, including all 48 Newton cases
%   reproduce('figures')  All seven manuscript figures
%   reproduce('full')     Complete tables and figures
%   reproduce('certify')  Tables plus all prose-claim and safeguard checks
%
% The default is 'quick'.  All CSV files are written to matlab/results.
% Figures go to tex/figures in the manuscript tree and to figures in a
% standalone code release.  Set CAVITATION_FIGURE_DIR to override this.

if(nargin<1 || isempty(mode)), mode='quick'; end
mode=lower(string(mode));
here=fileparts(mfilename('fullpath'));addpath(here);
outdir=fullfile(here,'results');

if(any(mode==["quick","tables","full","certify"]))
    fprintf('\n=== Tables 1--4: Reynolds verification ===\n')
    reproduce_reynolds_tables(outdir);
    fprintf('\n=== Tables 5--6: Stokes and model comparison ===\n')
    reproduce_stokes_core_tables(outdir);
end

if(any(mode==["tables","full","certify"]))
    fprintf('\n=== Table 7: constitutive laws ===\n')
    reproduce_constitutive_table(outdir);
    fprintf('\n=== Table 8: discretizations ===\n')
    reproduce_discretization_table(outdir);
    fprintf('\n=== Table 9: three-dimensional verification ===\n')
    reproduce_stokes_3d_table(outdir);
    fprintf('\n=== Table 10: smoothed Newton ===\n')
    reproduce_smoothed_newton_table(outdir);
end

if(mode=="certify")
    fprintf('\n=== Quantitative claims stated in prose ===\n')
    verify_manuscript_claims(outdir,'full');
    fprintf('\n=== Safeguarded Newton matrix ===\n')
    verify_newton_safeguard(outdir);
end

if(any(mode==["figures","full"]))
    fprintf('\n=== Reynolds pit fields ===\n')
    makepit
    makepitlambda
    makepitline
    fprintf('\n=== Stokes channel ===\n')
    makechannel
    fprintf('\n=== Reynolds--Stokes pit comparison ===\n')
    makepitcompare
    makepitfield
    fprintf('\n=== Three-dimensional Stokes verification ===\n')
    makestokes3d
end

if(~any(mode==["quick","tables","figures","full","certify"]))
    error('reproduce:mode','Unknown mode "%s"',mode)
end
fprintf('\nReproduction mode "%s" completed.\n',mode)
end
