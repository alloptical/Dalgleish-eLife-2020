function [o] = cat_structs(varargin)
% Concatenate one or more structures (with different fields) into a single
% structure

o = struct;
fields = cell(numel(varargin));
for v = 1:numel(varargin)
    fields{v} = fieldnames(varargin{v});
    for f = 1:numel(fields{v})
        o.(fields{v}{f}) = varargin{v}.(fields{v}{f});
    end
end

end

