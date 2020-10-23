function [out,flag] = pyb_return_trials(results,stim,var,varargin)
% quickly return trials of stim type stim and stim variation(s) var (i.e.
% can pass multiple variations for a given stim type and will return them
% all)
%
% if no additional arguments this will return a structure (with all fields
% intact). If you pass a field name as an additional argument this will
% return an array with just the values from that field.

flag = false(numel(results),1);
if ~isempty(var)
    for s = 1:numel(stim)
        for v = 1:numel(var)
            flag = flag | ([results(:).stim_type]'==stim(s) & [results(:).stim_var]'==var(v));
        end
    end
else
    for s = 1:numel(stim)
        flag = flag | [results(:).stim_type]'==stim(s);
    end
end
if ~isempty(varargin)
    if strcmp(varargin{1},'flag')
        out = flag;
    else
        out = [results(flag).(varargin{1})]';
    end
else
    out = results(flag);
end


end

