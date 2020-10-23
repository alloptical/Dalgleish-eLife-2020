function [out_dir] = backdir(in_dir,varargin)
% shifts a directory string back a directory (default) or by some number of
% directories (additional argument)

% Input = directory string to shift back
% Optional input: number of backwards steps to make

num = 1;
if ~isempty(varargin)
    num = varargin{1};
end
dirchunks = strsplit(in_dir,filesep);
dirchunks = dirchunks(cellfun(@(x) ~isempty(x),dirchunks));
out_dir = [filesep strjoin(dirchunks,filesep)];
for i = 1:num
    [out_dir,~] = fileparts(out_dir);
end
out_dir = [out_dir filesep];

end

