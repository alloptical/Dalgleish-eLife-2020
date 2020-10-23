%% Figure 2
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

scrsz = get(0,'Screensize');

% c: Target heatplots
figure('Position',[scrsz(3).*[0.12 0.5] scrsz(3:end).*[0.6 0.25]])
s = 7;
cLims = [-4 4];
toPlot = [5 50 200];

anyTargeted = any(cell2mat(targeted_flag(s,num_cells>0)),2);
nHeatmaps = numel(toPlot);
targHeatmaps = cell(1,numel(toPlot));
stimEpoch = stim_frame.ps:post_window_dff.ps(1)-1;
numCells = 200;
cellBar = round(0.1*numCells);

for i = 1:nHeatmaps
    targHeatmaps{i} = squeeze(mean(tw_traces_zs{s,num_cells == toPlot(i)}(anyTargeted,:,:),2));
    [~,order] = sort(mean(tw_response_zs{s,num_cells == toPlot(i)}(anyTargeted,:,:),2),'Descend');
    targHeatmaps{i} = targHeatmaps{i}(order,:);
    targHeatmaps{i}(:,stimEpoch) = 0;  
end
for i = 1:nHeatmaps
    subplot(1,nHeatmaps,i)
    imagesc(targHeatmaps{i}(1:numCells,:),cLims)
    set(gca,'XTick',[],'YTick',[])
    hold on
    xl = get(gca,'XLim'); yl = get(gca,'YLim');
    [xx,yy] = xy2vert(stimEpoch([1 end]) + [-0.5 0.5],yl);
    patch(xx,yy,'r','EdgeColor','none','FaceColor',[1 0.6 0.6])
    line(min(xl)+[0 imaging_rate_rounded],min(yl)*[1 1] + 0.5*range(yl),'Color',[0 0 0])
    line(min(xl)+[0 0],[min(yl)*[1 1] + 0.5*range(yl)] + [0 cellBar],'Color',[0 0 0])
    text(min(xl),[min(yl) + 0.51*range(yl)] + cellBar,num2str(cellBar))
end
colormap('redblue')
cb = colorbar('Location','southoutside','Orientation','horizontal','Ticks',round(sort([cLims 0]),2));
cb.Position = cb.Position .* [1 0.4 0.5 0.75];

% Psychometric curves
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.6 0.5]])

subplot(1,4,1)
mxSz = 500;
plotOptions = struct;
plotOptions.plotData = false;
plotOptions.CIthresh = false;
plotOptions.lineColor = [0.85 0.85 0.85];
plotOptions.fontSize = 10;
plotOptions.labelSize = 18;
plotOptions.extrapolLength = 5;
plotOptions.plotAsymptote = false;
plotOptions.plotThresh = false;
mkSz = normalise(nTrials)+0.01;
curvWid = 1+normalise(sum(nTrials,2));
curvCol = 0.75 + ((1 - normalise(sum(nTrials,2))) * 0.15);
thresh = zeros(num_sessions,1);
scatter(nCellsPsych(:),pResponse(:),mkSz(:)*mxSz,'k.','MarkerEdgeColor',[0.7 0.7 0.7])

% d: individuals
for s = 1:num_sessions
    hold on
    plotOptions.lineWidth = curvWid(s);
    plotPsych(numCellsBhvFitInd(s),plotOptions);
    thresh(s) = getThreshold(numCellsBhvFitInd(s),0.5);
    
    % formatting
    axis square
    set(gca,'XTick',[0.1 1 10 100],'XTickLabels',[0.1 1 10 100],'YTick',[0:0.25:1])
    ylabel('P(Lick)')
    xlim([min(nCellsPsych(:)) max(nCellsPsych(:))])
    
end

% d: group average
plotOptions.CIthresh = true;
plotOptions.dataColor = [0.4 0.4 0.4];
plotOptions.lineColor = [0.8 0 0];
plotOptions.lineWidth = 1.5;
plotOptions.fontSize = 12;
plotOptions.labelSize = 18;
plotOptions.extrapolLength = 0.05;
plotOptions.plotAsymptote = false;
plotOptions.plotThresh = false;
plotOptions.CIthresh = false;

plotPsych(numCellsBhvFit,plotOptions);
hold on
xlim([min(nCellsPsych(:)) max(nCellsPsych(:))])
set(gca,'XTick',[0.1 1 10 100],'XTickLabels',[0.1 1 10 100],'YTick',[0:0.25:1])
xlabel({'No. activated target neurons' '(log)'})
ylabel('P(Lick)')
title(['R2 = ' num2str(round(numCellsBhvFitR2,2))])
thisCI = abs(numCellsBhvFitPoints'-numCellsBhvFitPointsCI);
yVals = psignifit_x2y(numCellsBhvFit,log(pts2analyse));
miny = min(get(gca,'YLim'));
line(numCellsBhvFitPoints.*ones(2,3),[miny*ones(1,3) ; yVals],'Color',[0.8 0 0],'LineWidth',1.5,'LineStyle','-')
errorbar(numCellsBhvFitPoints,yVals,[],[],thisCI(:,1)',thisCI(:,1)',...
    'ko','Color',[0.8 0 0],'MarkerFaceColor',[0.8 0 0],'MarkerSize',5,'CapSize',0,'LineWidth',1.5)
xlim([1 max(get(gca,'Xlim'))])

ax = subplot(1,4,2);
hold on
h = plotSpread({numCellsBhvFitIndR2},'distributionMarkers','.','distributionColors',[0.6 0.6 0.6],'xValues',1.25);
boxplot(numCellsBhvFitIndR2,'Positions',0.75,'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
boxplotEdit('Median','Width',0.3,'Median','Color',[0.7 0.7 0.7],'Median','LineWidth',2,...
    'Box','Color',[0 0 0],'Box','LineWidth',10,'Whiskers','LineWidth',2,'Whisker','Color',[0.7 0.7 0.7])
xlim([0 2])
ylim([0.59 1.01])
axis square
ylabel('R2')
box off
set(gca,'XTick',[],'YTick',[0:0.2:1],'XColor','none','YAxisLocation','right')
ax.Position = ax.Position.*[1 1 0.4 0.8];

% e
subplot(1,4,3)
hold on
tmpx = repmat([1 2 3],num_sessions,1);
tmpy = numCellsBhvFitIndPoints;
plot(tmpx',tmpy','Color',[0.9 0.9 0.9],'Marker','.','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerSize',20)
boxplot(numCellsBhvFitIndPoints,'Positions',[1 2 3]+0.2,'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
boxplotEdit('Median','Width',0.2,'Median','Color',[0 0 0],'Median','LineWidth',2,...
    'Box','Color',[0 0 0],'Box','LineWidth',8,'Whiskers','LineWidth',2,'Whisker','Color',[0 0 0])
thisCI = abs(numCellsBhvFitPoints'-numCellsBhvFitPointsCI);
errorbar([1 2 3]+0.4,numCellsBhvFitPoints,thisCI(:,1),thisCI(:,2),'k.','Color',[0.8 0 0],'MarkerSize',20,'CapSize',0,'LineWidth',2)
axis_tight
axis square
box off
set(gca,'XTickLabels',[10 50 90])
xlabel('% point of curve')
ylabel('No. neurons')

% f
subplot(1,4,4)
hold on
h = plotSpread(numCellsBhvIndSlopes,'distributionColor',[0.6 0.6 0.6],'distributionMarkers','.');
set(h{1},'MarkerSize',15)
line([1.8 2.2],numCellsBhvFitSlope*[1 1],'Color',[1 0 0])
boxplot(numCellsBhvIndSlopes,'Positions',1.5,'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
boxplotEdit('Median','Width',0.3,'Median','Color',[0 0 0],'Median','LineWidth',2,...
    'Box','Color',[0 0 0],'Box','LineWidth',8,'Whiskers','LineWidth',2,'Whisker','Color',[0 0 0],'Outliers','YData',nan)
axis square
axis_tight
xlim([0.5 4])
box off
ylabel({'P(Lick) / neuron' 'at 50% point'})
set(gca,'XColor','none')
ylim([0 0.15])
set(gca,'YTick',[0:0.05:0.15])

%% Figure 2 - figure supplement 1h - j
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

responded = e_t_n;
respondedCatch = e_t_n_catch;

figure
subplot(2,2,1)
hold on
errorbar(mean(targeted,1),nanmean(allTargets,1),sem(allTargets,1),'k.-','Color',[0 0 0],'MarkerSize',20,'CapSize',0)
axis square
xlabel('No. target sites (log)')
ylabel('No. neurons in sites')
set(gca,'XScale','log','XTick',[1 10 25 50 100 200])
axis_tight

subplot(2,2,2)
hold on
errorbar(mean(targeted,1),nanmean(responded,1),sem(responded,1),'k.-','Color',[0 0 0],'MarkerSize',20,'CapSize',0)
errorbar(mean(targeted,1),nanmean(respondedCatch,1),sem(respondedCatch,1),'k.--','Color',[0.6 0.6 0.6],'MarkerSize',20,'CapSize',0)
axis square
xlabel('No. target sites (log)')
ylabel('No. responsive neurons in sites')
set(gca,'XScale','log','XTick',[1 10 25 50 100 200])
axis_tight

subplot(2,2,3)
hold on
x = targeted;
y = responded./allTargets;
yCatch = respondedCatch./allTargets;
line([min(targeted(targeted>0)) max(targeted(:))],0.05*[1 1])
errorbar(nanmean(targeted,1),nanmean(y,1),sem(y,1),'k.-','CapSize',0,'Color',[0 0 0],'MarkerSize',20)
errorbar(nanmean(targeted,1),nanmean(yCatch,1),sem(yCatch,1),'k.--','CapSize',0,'Color',[0.6 0.6 0.6],'MarkerSize',20)
axis_tight
axis square
xlabel('No. target sites (log)')
ylabel('N(responsive) / N(neurons in sites)')
set(gca,'XScale','log','XTick',[1 10 25 50 100 200])
axis_tight

subplot(2,2,4)
y = responded./allTargets;
yCatch = respondedCatch./allTargets;
categorical_plot([mean(y(:,~ng_var),2) mean(yCatch(:,~ng_var),2)],'x',{'Stim' 'Catch'})
ylabel('N(responsive) / N(neurons in sites)')

%% Figure 2 - figure supplement 2
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

nCells = e_t_n(:,~ng_var);
nCells(nCells<1) = nan;
nCells = log(nCells);
allRxnMu = cellfun(@(x) nanmean(x),tw_rxn(:,~ng_var));
allRxnSd = cellfun(@(x) nanstd(x),tw_rxn(:,~ng_var));
nRxn = cellfun(@(x) sum(~isnan(x)),tw_rxn(:,~ng_var));
flag = isnan(nCells) | isinf(nCells) | nRxn<3;

figure
subplot(1,2,1)
numCellsRxnMuFit = lin_fit_plot(nCells(~flag),allRxnMu(~flag),'EdgeColor','none');
ylim([-0.1 1.1])
axis_tight
ylabel({'Reaction time' 'mean (s)'})
xlabel({'No. activated target neurons' '(log)'})
set(gca,'XTick',log([1 10 100]),'XTickLabels',[1 10 100])

subplot(1,2,2)
numCellsRxnSdFit = lin_fit_plot(nCells(~flag),allRxnSd(~flag),'EdgeColor','none');
axis_tight
xlabel({'No. activated target neurons' '(log)'})
ylabel({'Reaction time' 's.d. (s)'})
set(gca,'XTick',log([1 10 100]),'XTickLabels',[1 10 100])

