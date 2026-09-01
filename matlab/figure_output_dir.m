function outdir=figure_output_dir(here)
% Select a figure directory for both the manuscript tree and code release.
override=getenv('CAVITATION_FIGURE_DIR');
if(~isempty(override))
    outdir=override;
elseif(isfolder(fullfile(here,'..','tex')))
    outdir=fullfile(here,'..','tex','figures');
else
    outdir=fullfile(here,'..','figures');
end
if(~exist(outdir,'dir')), mkdir(outdir); end
end
