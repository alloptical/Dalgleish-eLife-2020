function [all_bhv,animal_ids] = pyb_load_bhv_dir(bhv_dir,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

params.response_window     = 1.15;
params.sated_thold         = 0.85;
params.min_dur             = 0.5;
params.min_lick_t          = 0.3;
params.remove_autoreward   = false;
params.remove_ambiguous    = true;
if ~isempty(varargin)
    tmp = varargin{1};
    f = fieldnames(tmp);
    for i = 1:numel(f)
        params.(f{i}) = tmp.(f{i});
    end
end

[animal_dirs,animal_ids] = return_animal_dirs(bhv_dir);
all_bhv = [];
n = 0;
for i = 1:numel(animal_dirs)
    paths = return_bhv_paths(animal_dirs{i});
    if ~isempty(paths)
        n = n+1;
        all_bhv(n).animal_id = animal_ids{i};
        all_bhv(n).bhv_paths = paths;
        num_sessions = numel(all_bhv(n).bhv_paths);
        [all_bhv(n).params,all_bhv(n).stim_vars,all_bhv(n).bhv_init,all_bhv(n).bhv] = deal(cell(num_sessions,1));
        for j = 1:num_sessions
            load(all_bhv(n).bhv_paths{j},'results')
            all_bhv(n).params{j} = struct;
            
            [all_bhv(n).bhv_init{j},all_bhv(n).stim_vars{j},p] = ...
                pyb_parse_results(results,'response_window',params.response_window,'remove_autoreward',params.remove_autoreward);
            
            all_bhv(n).params{j} = cat_structs(all_bhv(n).params{j},p);
            
            [all_bhv(n).bhv{j},p] = ...
                proc_bhv_lr(all_bhv(n).bhv_init{j},...
                'sated_thold',params.sated_thold,'min_dur',params.min_dur,...
                'min_lick_t',params.min_lick_t,'response_window',params.response_window,...
                'remove_ambiguous',params.remove_ambiguous);
            
            all_bhv(n).params{j} = cat_structs(all_bhv(n).params{j},p);
        end
        try
            all_bhv(n).bhv_cat = cell2mat(all_bhv(n).bhv);
        catch
            all_bhv(n).bhv_cat = [];
        end
    end
end

end

