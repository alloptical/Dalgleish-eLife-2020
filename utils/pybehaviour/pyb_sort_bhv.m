function [bhv,order,stimvar_trials,stim_vars] = pyb_sort_bhv(bhv,varargin)
% sort pybehaviour structure, which is in order of trials delivered, but
% order of stimulus types and variations defined by a stimvars matrix (row
% 1 = stimulus types, row 2 = stimulus variations)

stim_vars = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'StimVars')
        stim_vars = varargin{v+1};
    end
end
if isempty(stim_vars)
    stim_vars = pyb_return_stimvars(bhv);
end

trial_order = [[bhv(:).stim_type]' [bhv(:).stim_var]'];
order = [];
stimvar_trials = cell(1,size(stim_vars,2));
for s = 1:size(stim_vars,2)
    order = [order ; find(min(trial_order==stim_vars(:,s)',[],2)==1)];
    stimvar_trials{s} = find(min(trial_order==stim_vars(:,s)',[],2)==1);
end
bhv = bhv(order);
end

