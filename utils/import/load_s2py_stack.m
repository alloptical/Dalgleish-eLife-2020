function [f,np,centroids,stat,ops,spks] = load_s2py_stack(s2p_path)

%% load_s2py_stack

% Henry Dalgleish 20191023
% Function to flexibly import Python Suite2P .mat files into matlab.

% Note this function will automatically filter by iscell

% Inputs: --> single file path: either a "combined" file (multiple planes
% in one Fall.mat) or a single plane experiment --> cell array: each cell
% corresponds to a plane from a single experiment containing either an
% Fall.mat file path or a Suite2P "planeN" directory containing a Fall.mat
% file.

% Outputs: self-explanatory except --> with single file path input the ops
% variable corresponds EITHER to that plane, or the combined Fall.mat data
% --> with n-dimensional cell array input ops variable is a structure of
% dimensions n where each entry is for a given plane

%% main

% for "combined" file
if ischar(s2p_path)
    [f,np,centroids,stat,ops,spks] = load_s2p_file(s2p_path);
    
% for individual files
elseif iscell(s2p_path)
    f = [];
    np = [];
    centroids = [];
    stat = [];
    iplanes = [];
    for i = 1:numel(s2p_path)
        if exist(s2p_path{1},'file')==2
            [temp_f,temp_np,temp_centroids,temp_stat,temp_ops] = load_s2p_file(s2p_path{i});
        elseif exist(s2p_path{1},'dir') == 7
            d = dir(s2p_path{i});
            d = {d(:).name};
            in_file = fullfile(s2p_path{i},d{strcmp(d,'Fall.mat')});
            [temp_f,temp_np,temp_centroids,temp_stat,temp_ops] = load_s2p_file(in_file);
        end
        f = [f ; temp_f];
        np = [np ; temp_np];
        temp_centroids(:,end) = temp_centroids(:,end)*i;
        centroids = [centroids ; temp_centroids];
        stat = [stat temp_stat];
        ops(i) = temp_ops;
        iplanes = [iplanes ; temp_centroids(:,end)];
    end
    for i = 1:numel(stat)
        stat(i).iplane = iplanes(i);
    end
    
else
    fprintf('>>>>>>>>>>>>>>>>>> Warning: Input not recognised <<<<<<<<<<<<<<<<<<<<\nPass either a single combined Fall.mat file or a cell of directories\nor files, one for each plane\n')
end

end

%% subfunctions

function [f,np,centroids,stat,ops,spks] = load_s2p_file(s2p_path)
iscell = [];
load(s2p_path);
ic = iscell(:,1) == 1;
stat = stat(ic );
stat = [stat{:}];
f = F(ic,:);
np = Fneu(ic,:);
spks = spks(ic,:);

% for suite2py combined outputs need to reindex x and y (they are currently
% in the "stitched space" used in the suite2py gui).
if isfield(stat,'xpix_original')
    for i = 1:numel(stat)
        tmp = [stat(i).xpix ; stat(i).ypix];
        stat(i).xpix = stat(i).xpix_original;
        stat(i).ypix = stat(i).ypix_original;
        stat(i).xpix_offset = tmp(1,:);
        stat(i).ypix_offset = tmp(2,:);
    end
else
    for i = 1:numel(stat)
        [x,y] = ind2sub(size(ops.refImg),stat(i).ipix);
        stat(i).xpix_offset = stat(i).xpix;
        stat(i).ypix_offset = stat(i).ypix;
        stat(i).xpix = x;
        stat(i).ypix = y;
    end
end
% for individual plane files add iplane field and set to 0
if ~isfield(stat,'iplane')
    for i = 1:numel(stat)
        stat(i).iplane = 0;
    end
end
centroids = [round(cellfun(@(x) median(x),[{stat(:).xpix}' {stat(:).ypix}'])) [stat.iplane]'+1];

% loop through each plane creating roi ims/hsvs
planes = unique([stat(:).iplane]);
n_planes = numel(planes);
ops.roi_ims = cell(n_planes,1);
try
    lims_yx = [ops.dy ops.dx];
catch
    lims_yx = [ops.Ly ops.Lx];
end
for p = 1:n_planes
    these_cells = centroids(:,3)==p;
    ops.roi_ims{p} = zeros(lims_yx(1),lims_yx(2));
    ipix = sub2ind(lims_yx,[stat(these_cells).ypix],[stat(these_cells).xpix]);
    ops.roi_ims{p}(ipix) = [stat(these_cells).lam];
    ops.roi_hsv{p} = ...
        s2p_rois2image({stat(these_cells).xpix},{stat(these_cells).ypix},{stat(these_cells).lam},lims_yx);
end
end

