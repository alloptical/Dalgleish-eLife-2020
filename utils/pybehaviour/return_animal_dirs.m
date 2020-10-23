function [paths,animal_ids,d] = return_animal_dirs(path)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

d = dir(path);
filter = '[a-zA-Z]{1,2}[0-9]{3}$';
d = d(cellfun(@(x) ~isempty(regexp(x,filter,'ONCE')),{d(:).name}));
[paths,animal_ids] = deal(cell(numel(d),1));
for i = 1:numel(d)
    paths{i} = fullfile(d(i).folder,d(i).name);
    animal_ids{i} = d(i).name;
end

end

