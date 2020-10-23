%% Setup
% Ensures repository is added to matlab path 
%
% NB by default these figure scripts will load the intermediate processed
% data included in this repo, however this can be rederived from raw data
% downloaded from: 
% https://doi.org/10.6084/m9.figshare.13128950
%
% see DalgleishHausser2020_analysis_main.m for more details
%
% Note this plotting script must be in the DalgleishHausser2020 base
% directory

[repo_path,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
run([repo_path filesep 'DalgleishHausser2020_setup.m'])

%% Figure 1
run('DalgleishHausser2020_Fig1_plots')

%% Figure 2
run('DalgleishHausser2020_Fig2_plots')

%% Figure 3
% Figure 3 figure supplement 1i-n require raw Figshare data (see above). If
% you have this and want to plot these panels then set the import_raw to
% true. If false it'll skip these panels.
import_raw = false; 
run('DalgleishHausser2020_Fig3_plots')

%% Figure 4
run('DalgleishHausser2020_Fig4_plots')
