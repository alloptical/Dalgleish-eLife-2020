%% Setup
% Ensures repository is added to matlab path
% Note this analysis script must be in the DalgleishHausser2020 base
% directory

[repo_path,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
run([repo_path filesep 'DalgleishHausser2020_setup.m'])

%% Raw imaging data import and processing
% Note this requires data downloaded from:
% https://doi.org/10.6084/m9.figshare.13128950
%
% unzip this data using unzipDirectories.m in this repo and in the below
% script set:
% base_dir = 'path/to/data'
% where this points to a directory of unzipped animal session directories:
%
% --- ~/DalgleishHausser2020Repo/rawDataDirectory/
%     |--- 20191203_L541/
%     |    --- experimentDataFile1
%     |    --- experimentDataFile2
%     |--- 20191009_L544/
%     |    --- experimentDataFile1
%     |    --- experimentDataFile2
%     . etc.

run('DalgleishHausser2020_importProcessing')

% Cells below this can be run either with the data downloaded and imported
% above, or with saved versions of intermediate processing steps included
% in this github repository:
% https://github.com/alloptical/Dalgleish-eLife-2020
%
% (note you can skip this cell if you don't have raw data)

%% Threshold estimation procedure
run('DalgleishHausser2020_threshold_analysis')

%% Analyse trial-wise data for Figures 2 - 3
run('DalgleishHausser2020_Fig2_3_4_analysis')

%% Analyse 1P and 2P training data
run('DalgleishHausser2020_Fig1_analysis')
