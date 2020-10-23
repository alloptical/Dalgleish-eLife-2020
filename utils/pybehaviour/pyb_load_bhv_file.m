function [bhv] = pyb_load_bhv_file(bhv_file,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

params.response_window     = 2;
params.sated_thold         = 0.7;
params.min_dur             = 0.5;
params.min_lick_t          = 0.15;
params.remove_autoreward   = false;
params.remove_ambiguous    = true;
if ~isempty(varargin)
    tmp = varargin{1};
    f = fieldnames(tmp);
    for i = 1:numel(f)
        params.(f{i}) = tmp.(f{i});
    end
end

load(bhv_file)
[bhv.bhv_init,bhv.stim_vars,p] = pyb_parse_results(results,'response_window',params.response_window,'remove_autoreward',params.remove_autoreward);
bhv.params = struct;
bhv.params = cat_structs(bhv.params,p);
[bhv.bhv_proc,p] = proc_bhv_lr(bhv.bhv_init,'sated_thold',params.sated_thold,'min_dur',params.min_dur,'min_lick_t',params.min_lick_t,'response_window',params.response_window,'remove_ambiguous',params.remove_ambiguous);
bhv.params = cat_structs(bhv.params,p);

end

