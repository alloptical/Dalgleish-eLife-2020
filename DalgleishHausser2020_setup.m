% Ensure this script is in the repository's base directory, i.e.:
% ~/DalgleishHausser2020Repo/DalgleishHausser2020_setup.m

[repo_path,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(repo_path)
addpath(genpath(repo_path))

if ~exist('psignifit-HD')
    warning('Repository path not set up correctly. Ensure DalgleishHausser2020_setup.m is in the parent directory and rerun:\n%s',...
            '~/DalgleishHausser2020Repo/DalgleishHausser2020_setup.m')
else
    fprintf(['\nRepository paths added successfully!\n\n'])
end