function stim_var_targets = load_tpbs_targets(path)

% Loads photostimulation Points files used in behavioural experiments and
% syncs them with the behavioural VarFile to return xyz positions of
% photostimulation sites on each trial type.
%
% Input = path to directory containing Points and VarFile matlab files

% parse data
tpbs_path = return_fullfile(path,'_VarFile_');
targets_path = sort_nat(return_fullfile(path,'_Points.mat'));

% load data
load(tpbs_path);
num_targets = numel(targets_path);
naparms = cell(num_targets,1);
for i = 1:num_targets
    naparms{i} = load(targets_path{i});
end

% sort phasemasks into protocol rows
var_phasemasks = vf.phasemasks.phasemasks_details(:,1);
[targets_xyz,naparm_str] = deal(cell(num_targets,1));
[group_id,naparm_id] = deal(zeros(num_targets,1));
for i = 1:numel(var_phasemasks)
    idx = strfind(var_phasemasks{i},'NAPARM_');
    naparm_str{i} = var_phasemasks{i}(idx:idx+9);
    naparm_id(i) = find(contains(targets_path,naparm_str{i}));
    group_id(i) = str2double(var_phasemasks{i}(1:3));
    
    these_points = naparms{naparm_id(i)}.points.Group==group_id(i);
    targets_xyz{i} = [naparms{naparm_id(i)}.points.X(these_points)' ...
                      naparms{naparm_id(i)}.points.Y(these_points)' ...
                      naparms{naparm_id(i)}.points.Z(these_points)'];
end

% for each stim/var find the appropriate targets according to protocol
tbl = vf.protocol.protocol_table(:,3);
num_rows = size(tbl,1);
phasemask_ids = cell(num_rows,1);
stim_var_targets = struct;
for i = 1:num_rows
    stim_var_targets(i).stim_type = vf.protocol.stim_var_mat(i,1);
    stim_var_targets(i).stim_var = vf.protocol.stim_var_mat(i,2);
    if isnan(tbl{i})
        phasemask_ids{i} = [];
    else
        phasemask_ids{i} = cellfun(@(x) str2double(x),strsplit(tbl{i},' '));
        stim_var_targets(i).targets_xyz = [];
        for j = 1:numel(phasemask_ids{i})
            stim_var_targets(i).targets_xyz = [stim_var_targets(i).targets_xyz ; targets_xyz{phasemask_ids{i}(j)}];
        end
    end
end


end

