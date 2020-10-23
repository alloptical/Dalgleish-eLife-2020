function [stim_vars] = pyb_return_stimvars(bhv)
% return stimulus types and variations for a given session bhv

stim_types = double(unique([bhv(:).stim_type]));
stim_vars = [];
for i = 1:numel(stim_types)
    vars = double(unique(pyb_return_trials(bhv,stim_types(i),[],'stim_var')));
    stim_vars = [stim_vars [stim_types(i) * ones(1,numel(vars)) ; vars']];
end

end

