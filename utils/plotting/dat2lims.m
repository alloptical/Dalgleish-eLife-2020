function [lims] = dat2lims(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[mi,ma] = deal(nan);
factor = 0.1;
skip = [];
for v = 1:numel(varargin)
    if ~ismember(v,skip)
        if strcmpi(varargin{v},'factor')
            factor = varargin{v+1};
            skip = [skip v+1];
        else
            mi = min([mi ; varargin{v}(:)]);
            ma = max([ma ; varargin{v}(:)]);
        end
    end
end
multiplier = [1 1] + factor(:)';
lims = [mi ma] .* [multiplier(1) multiplier(2)];

end

