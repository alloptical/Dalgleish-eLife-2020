function [stim_vars] = pyb_lick_raster(bhv,varargin)
% Create a lick raster from imported/processed pybehaviour data.
% To sort rasters by specific order of stimulus types and variations pass a
% stim var matrix where first row is stim types and second row is variations
%
% bhv is a pybehaviour session structure returned by pyb_load_bhv_file.m
% and associated functions


%% Setup
stim_vars = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'StimVars')
        stim_vars = varargin{v+1};
    end
end

%%
% Figure
figure('Position',[708   512   150   515])
rw_end = bhv(1).parameters.responseWindow;
if isfield(bhv(1).parameters,'params')
    rw_start = bhv(1).parameters.params.min_lick_t;
else
    rw_start = 0;
end

offset = 15;
if isempty(stim_vars)
    %sorted_bhv = pyb_sort_bhv(bhv);
    [licks_x,licks_y,rxn_x,rxn_y,c,trial_type] = pyb_return_raster(bhv,'sort',true);
else
    %sorted_bhv = pyb_sort_bhv(bhv,'StimVars',stim_vars);
    [licks_x,licks_y,rxn_x,rxn_y,c,trial_type] = pyb_return_raster(bhv,'sort',true,'StimVars',stim_vars);
end

%cumulative_trials = [0 ; cellfun(@(x) numel(x),trial_type)];
start_stop = zeros(numel(trial_type),2);
for i = 1:numel(trial_type)
    o = (i-1)*offset;
    licks_y(trial_type{i}) = cellfun(@(x) x+o,licks_y(trial_type{i}),'UniformOutput',0);
    rxn_y(trial_type{i}) = rxn_y(trial_type{i}) + o;
    start_stop(i,:) = [trial_type{i}(1) trial_type{i}(end)] + o;
end
xl = [0 rw_end+1.2];
yl = [1 max(start_stop(:))];
hold on
rw_patch_x = repmat([rw_start rw_start rw_end rw_end]',1,size(start_stop,1));
rw_patch_y = [start_stop fliplr(start_stop)]';
patch(rw_patch_x,rw_patch_y,'r','FaceColor',[0.9 0.9 0.9],'EdgeColor','none','LineWidth',2)
plot(xl,start_stop(:)*[1 1],'Color',[0.6 0.6 0.6],'LineWidth',2)

lx = cell2mat(licks_x);
ly = cell2mat(licks_y);
flag = lx>rw_end+1;
lx(flag) = [];
ly(flag) = [];
flag = rxn_x>rw_end+1;
rxn_x(flag) = [];
rxn_y(flag) = [];
scatter(lx,ly,30,'ko','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7])
scatter(rxn_x,rxn_y,60,'ko','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.2 0.2 0.2])

plot([0 0],start_stop,'Color',[0 0 0],'LineWidth',5)
plot(rw_end*[1 1],start_stop,'Color',[0 0 0],'LineWidth',2)

for i = 1:numel(trial_type)
    imagesc(2.55,start_stop(i,1),c(trial_type{i}),[0 1])
    patch([2.05 2.05 2.2 2.2],[start_stop(i,:) fliplr(start_stop(i,:))],'r','FaceColor','none','EdgeColor',[0 0 0],'LineWidth',2)
end
colormap('gray')

xlim(xl)
ylim(yl)
set(gca,'XTick',[0:1:xl(end)],'YTick',start_stop(:,end),'YTickLabels',diff(start_stop,[],2)+1,'Layer','top')
xlabel('Reaction time (s)')


end


