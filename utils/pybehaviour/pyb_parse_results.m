function [bhv,stim_vars,params] = pyb_parse_results(results,varargin)
% Lower level import function - converts raw pybehaviour results cell to
% bhv struct

params.response_window = [];
params.remove_autoreward = false;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'remove_autoreward')
        params.remove_autoreward = varargin{v+1};
    elseif strcmpi(varargin{v},'response_window')
        params.response_window = varargin{v+1};
    end
end

results = results(cellfun(@(x) ~isempty(fieldnames(x)),results));
bhv = [results{:}]';
for t = 1:numel(bhv)
    responses = bhv(t).responses;
    if ~isempty(responses)
        responses = bhv(t).responses(1,:);
        if isempty(params.response_window)
            rw = bhv(t).parameters.responseWindow;
        else
            rw = params.response_window;
        end
        rxntime = responses(find(responses>0 & responses<=rw,1,'First'));
        if ~isempty(rxntime)
            bhv(t).rxntime = rxntime;
        else
            bhv(t).rxntime = nan;
        end
    else
        bhv(t).rxntime = nan;
    end
    bhv(t).responded = ~isnan(bhv(t).rxntime);
    bhv(t).istrial = true;
    bhv(t).trial_num = t;
end
if params.remove_autoreward
    bhv([bhv(:).auto_reward]==1) = [];
end

stim_types = double(unique([bhv(:).stim_type]));
stim_vars = [];
for i = 1:numel(stim_types)
    vars = double(unique(pyb_return_trials(bhv,stim_types(i),[],'stim_var')));
    stim_vars = [stim_vars [stim_types(i) * ones(1,numel(vars)) ; vars']];
end

end

