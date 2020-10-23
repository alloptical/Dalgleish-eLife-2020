function [licks_x,licks_y,rxn_x,rxn_y,c,trial_type] = pyb_return_raster(bhv,varargin)
% lower level behaviour raster function

stim_vars = [];
sort_flag = false;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'StimVars')
        stim_vars = varargin{v+1};
    elseif strcmpi(varargin{v},'sort')
        sort_flag = varargin{v+1};
    end
end

min_lick_t = bhv(1).parameters.params.min_lick_t;
if isempty(stim_vars)
    stim_vars = pyb_return_stimvars(bhv);
end
if sort_flag
    bhv = pyb_sort_bhv(bhv,'StimVars',stim_vars);
end

trial_order = [[bhv(:).stim_type]' [bhv(:).stim_var]'];
licks_x     = {bhv(:).responses};
responses   = [bhv(:).responded]'==1 & [bhv(:).rxntime]'>min_lick_t;
num_trials  = numel(responses);
rxn_x       = [bhv(:).rxntime]';
rxn_y       = 1:num_trials;

no_lick = cellfun(@(x) isempty(x),licks_x);
licks_x(~no_lick) = cellfun(@(x) x(1,x(1,:)>0),licks_x(~no_lick),'UniformOutput',0);
licks_y = cell(1,numel(licks_x));
for i = 1:numel(licks_y)
    licks_y{i} = i*ones(size(licks_x{i}));
end
licks_x(no_lick) = {nan};
licks_y(no_lick) = {nan};
c = ones(num_trials,1);
c(responses==1) = 0.2;
c(responses==0) = 0.9;

trial_type = cell(size(stim_vars,2),1);
for i = 1:size(stim_vars,2)
    trial_type{i} = find(min(trial_order == stim_vars(:,i)',[],2)==1);
end
    
end

