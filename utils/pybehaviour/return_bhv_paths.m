function [paths,d] = return_bhv_paths(path)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[~,animal_id] = fileparts(path);
d = dir([path filesep '*.mat']);
d = d(cellfun(@(x) ~isempty(regexp(x,[animal_id '_\d{8}_\d{6}.mat'],'ONCE')),{d(:).name}));
paths = cell(numel(d),1);
for i = 1:numel(d)
    paths{i} = fullfile(d(i).folder,d(i).name);
end

end

