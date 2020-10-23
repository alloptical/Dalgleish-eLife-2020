function [full_file] = return_fullfile(directory2search,varargin)
% Returns full filepath string for files in directory2search that match
% substring and/or regexp as a cell.
%
% Can combine raw string searches with regexp. Regexp is searched first,
% then raw strings 

d = dir(directory2search);
file_list = {d(:).name}';
flag = true(numel(file_list),1);

re_flags = find(cellfun(@(x) strcmpi(x,'re'),varargin));
re_varargin = varargin(re_flags+1);
varargin([re_flags(:) ; re_flags(:)+1]) = [];

% process regexp
for v = 1:numel(re_varargin)
    if ~iscell(re_varargin{v})
        flag = flag & cellfun(@(x) ~isempty(regexp(x,re_varargin{v},'ONCE')),file_list);
    elseif iscell(re_varargin{v})
        for f = 1:numel(re_varargin{v})
            flag = flag & cellfun(@(x) ~isempty(regexp(x,re_varargin{v}{f},'ONCE')),file_list);
        end
    end
end

% process raw string comparisons
for v = 1:numel(varargin)
    if ~iscell(varargin{v})
        flag = flag & cellfun(@(x) contains(x,varargin{v}),file_list);
    elseif iscell(varargin{v})
        for f = 1:numel(varargin{v})
            flag = flag & cellfun(@(x) contains(x,varargin{v}{f}),file_list);
        end
    end
end

if sum(flag) == 1
    full_file = [directory2search filesep file_list{flag}];
elseif sum(flag) > 1
    idcs = find(flag);
    full_file = cell(numel(idcs),1);
    for i = 1:numel(idcs)
        full_file{i} = [directory2search filesep file_list{idcs(i)}];
    end
else
    full_file = [];
end

end

