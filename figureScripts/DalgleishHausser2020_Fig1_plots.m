%% Figure 1
% stims are in stim_type/stim_var order, [5 6 7 7 ; 0 0 0 1]
% this corresponds to [1P Catch 2P_200 2P_100]
% these are re-ordered into [Catch 2P 1P] below for plotting

% load
load('DalgleishHausser2020_All2PLearning_analysed.mat')

scrsz = get(0,'Screensize');
% example behaviour session
pyb_lick_raster(bhv_eg.bhv_proc,'StimVars',[6 7 5 ; 0 0 0])
set(gcf,'Position',[0 0 scrsz(3:end).*[0.12 0.5]])

% group plots
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.12 0.5]])
%figure('Position',[508 424 281 531])

% 1e
subplot(3,1,1)
yPr = lastPr(:,[2 3 1]);
yPr = yPr(~any(isnan(yPr(:,1:3)),2),:);
plot(yPr','Color',[0.9 0.9 0.9],'Marker','.','MarkerSize',15)
hold on
errorbar(1:size(yPr,2),mean(yPr,1),sem(yPr,[],1),'k.','MarkerFaceColor',[0 0 0],'LineStyle','-','LineWidth',1.5,'CapSize',0,'MarkerSize',30)
box off
ylabel('P(Lick)')
set(gca,'XTick',1:size(yPr,2),'XTickLabels',{'Catch' '2P' '1P'},'YTick',[0:0.5:1])
axis_tight
axis square
xtickangle(30)
[yPr_p,~,yPr_stats] = friedman(yPr,1,'off');
yPr_mc = multcompare(yPr_stats,'CType','bonferroni','Display','off');

% 1f
subplot(3,1,2)
hold on
yMu = lastRxnMean(:,[2 3 1]);
yMu = yMu(~any(isnan(yMu),2),:);
plot(yMu','Color',[0.9 0.9 0.9],'Marker','.','MarkerSize',15)
errorbar(1:size(yMu,2),mean(yMu,1),sem(yMu,[],1),'k.','MarkerFaceColor',[0 0 0],'LineStyle','-','LineWidth',1.5,'CapSize',0,'MarkerSize',30)
box off
ylabel({'Mean reaction' 'time (s)'})
set(gca,'XColor','none')
axis_tight
set(gca,'XTick',1:size(yMu,2),'XTickLabels',{'Catch' '2P' '1P'})
ylim([0 1])
axis square
xtickangle(30)
% stats
names = {'Catch' '2P' '1P'};
group = reshape(repmat([1 2 3],size(yMu,1),1),[],1);
group = names(group);
t = table(group(:),yMu(:),'VariableNames',{'group' 'meas'});
yMu_rm = fitrm(t,'meas~group');
yMu_mc = multcompare(yMu_rm,{'group'},'ComparisonType','Bonferroni');

% 1g
subplot(3,1,3)
hold on
ySD = lastRxnSd(:,[2 3 1]);
ySD = ySD(~any(isnan(ySD),2),:);
plot(ySD','Color',[0.9 0.9 0.9],'Marker','.','MarkerSize',15)
errorbar(1:size(ySD,2),mean(ySD,1),sem(ySD,[],1),'k.','MarkerFaceColor',[0 0 0],'LineStyle','-','LineWidth',1.5,'CapSize',0,'MarkerSize',30)
box off
ylabel({'Reaction time' 'variability (s)'})
axis_tight
set(gca,'XTick',1:size(ySD,2),'XTickLabels',{'Catch' '2P' '1P'},'YTick',[0:0.2:0.4])
ylim([0 max(get(gca,'YLim'))])
axis square
xtickangle(30)
set(gca,'XColor','none')
% stats
names = {'Catch' '2P' '1P'};
group = reshape(repmat([1 2 3],size(ySD,1),1),[],1);
group = names(group);
t = table(group(:),ySD(:),'VariableNames',{'group' 'meas'});
ySD_rm = fitrm(t,'meas~group');
ySD_mc = multcompare(ySD_rm,{'group'});

%% Figure 1 - figure supplement 3
load('DalgleishHausser2020_All1PLearning_analysed.mat')

numTrials2Plot = 300;

% c
figure('Position',get(0,'Screensize').*[0 0 0.5 0.5])
subplot(3,3,1)
hold on
axis_tight
shadedErrorBar(1:numTrials2Plot,nanmean(slidingPr(:,1:numTrials2Plot),1),sem(slidingPr(:,1:numTrials2Plot),1),{'Color',[1 0.6 0]})
shadedErrorBar(1:numTrials2Plot,nanmean(slidingPrCatch(:,1:numTrials2Plot),1),sem(slidingPrCatch(:,1:numTrials2Plot),1),{'Color',[0.6 0.6 0.6]})
xlim([0-10 numTrials2Plot+10])
ylim([-0.1 1.1])
set(gca,'XTick',[0:100:numTrials2Plot])
ylabel({'P(Lick)' 'P(Lick) before auto-reward'})
axis square
xlabel('Trial')

% d
subplot(3,3,2)
hold on
line([0 numTrials2Plot],[0.5 0.5],'Color',[0.4 0.4 0.4])
shadedErrorBar(1:numTrials2Plot,nanmean(slidingRxnCatch(:,1:numTrials2Plot),1),sem(slidingRxnCatch(:,1:numTrials2Plot),1),{'Color',[0.6 0.6 0.6]})
shadedErrorBar(1:numTrials2Plot,nanmean(slidingRxn(:,1:numTrials2Plot),1),sem(slidingRxn(:,1:numTrials2Plot),1),{'Color',[1 0.6 0]})
axis_tight
xlim([0-10 numTrials2Plot+10])
set(gca,'XTick',[0:100:numTrials2Plot])
ylabel({'Reaction time mean' '(s)'})
axis square
xlabel('Trial')

% e
subplot(3,3,3)
hold on
shadedErrorBar(1:numTrials2Plot,nanmean(slidingRxnSdCatch(:,1:numTrials2Plot),1),sem(slidingRxnSdCatch(:,1:numTrials2Plot),1),{'Color',[0.6 0.6 0.6]})
shadedErrorBar(1:numTrials2Plot,nanmean(slidingRxnSd(:,1:numTrials2Plot),1),sem(slidingRxnSd(:,1:numTrials2Plot),1),{'Color',[1 0.6 0]})
axis_tight
xlim([0-10 numTrials2Plot+10])
set(gca,'XTick',[0:100:numTrials2Plot])
ylabel({'Reaction time s.d.' '(s)'})
axis square
xlabel('Trial')

% f
subplot(3,3,4)
col = copper(140);
col = col(end-101:end,:);
hold on
x = 1:size(sw_mean_dp_cat,2);
yyaxis left
errorbar(x,nanmean(sw_mean_dp_cat,1),sem(sw_mean_dp_cat,1),'k.-','CapSize',0,'MarkerSize',18)
axis_tight;
ylabel('Mean d-prime')
yyaxis right
yMu = nanmean(sw_mean_power_cat,1);
yErr = sem(sw_mean_power_cat,1);
plot(x,yMu,'Color',col(7*10,:))
for i = 1:numel(x)
    errorbar(x(i),yMu(i),yErr(i),'k.','MarkerSize',18,'CapSize',0,...
        'Color',col(ceil(yMu(i))*10,:))
end
axis_tight;
set(gca,'YColor',col(7*10,:))
ylabel('Mean power (mW)')
xlim([0 max(x)+1])
axis square
box off
set(gca,'XTick',x)
xlabel('Session')

% g
subplot(3,3,5)
col = copper(140);
col = col(end-101:end,:);
hold on
x = 1:size(sw_min_dp_cat,2);
yyaxis left
errorbar(x,nanmean(sw_min_dp_cat,1),sem(sw_min_dp_cat,1),'k.-','CapSize',0,'MarkerSize',18)
axis_tight;
ylabel('d-prime at min power')
yyaxis right
yMu = nanmean(sw_min_power_cat,1);
yErr = sem(sw_min_power_cat,1);
plot(x,yMu,'Color',col(7*10,:))
for i = 1:numel(x)
    errorbar(x(i),yMu(i),yErr(i),'k.','MarkerSize',18,'CapSize',0,...
        'Color',col(ceil(yMu(i))*10,:))
end
axis_tight;
set(gca,'YColor',col(7*10,:))
ylabel('Min power (mW)')
xlim([0 max(x)+1])
axis square
box off
set(gca,'XTick',x)
xlabel('Session')

% h
h = subplot(3,3,6);
c = [0.7 0.7 0.7 ; 0.5 0.5 0.5 ; 0 0 0];
for d = 1:maxSessionsHighPower 
    hold on
    errorbar(highPowers*1000,nanmean(crossDayHighPsychsPr(:,:,d),1),sem(crossDayHighPsychsPr(:,:,d),1),'k.-','MarkerSize',20,'Color',c(d,:),'CapSize',0)
end
axis_tight
axis square
set(gca,'XTick',sort(highPowers*1000),'YTick',[0:0.25:1])
xlabel('Power (µW)')
ylabel('P(Lick)')

% i
subplot(3,3,7);
c = [0.7 0.7 0.7 ; 0.5 0.5 0.5 ; 0 0 0];
for d = 1:3
    hold on
    errorbar(lowPowers*1000,nanmean(crossDayLowPsychsPr(:,:,d),1),sem(crossDayLowPsychsPr(:,:,d),1),'k.-','MarkerSize',20,'Color',c(d,:),'CapSize',0)
end
[xl,~] = axis_tight;
axis square
xlim(xl)
set(gca,'XTick',sort(lowPowers*1000),'YTick',[0:0.25:1])
xlabel('Power (µW)')
ylabel('P(Lick)')

% j
subplot(3,3,8);
c = [0.7 0.7 0.7 ; 0.5 0.5 0.5 ; 0 0 0];
for d = 1:3 
    hold on
    errorbar(lowPowers*1000,nanmean(crossDayLowPsychsRxnMean(:,:,d),1),sem(crossDayLowPsychsRxnMean(:,:,d),1),'k.-','MarkerSize',20,'Color',c(d,:),'CapSize',0)
end
axis_tight;
axis square
set(gca,'XTick',sort(lowPowers*1000),'YTick',[0:0.25:1],'YLim',[0 1])
xlabel('Power (µW)')
ylabel('Reaction time mean (s)')

% k
subplot(3,3,9);
c = [0.7 0.7 0.7 ; 0.5 0.5 0.5 ; 0 0 0];
for d = 1:3
    hold on
    errorbar(lowPowers*1000,nanmean(crossDayLowPsychsRxnSd(:,:,d),1),sem(crossDayLowPsychsRxnSd(:,:,d),1),'k.-','MarkerSize',20,'Color',c(d,:),'CapSize',0)
end
axis_tight;
axis square
set(gca,'XTick',sort(lowPowers*1000),'YTick',[0:0.25:1],'YLim',[0 1])
xlabel('Power (µW)')
ylabel('Reaction time s.d. (s)')

%% Figure 1 - figure supplement 5
load('DalgleishHausser2020_All2PLearning_analysed.mat')
learning_1p = load('DalgleishHausser2020_All1PLearning_analysed.mat','nSessionsStaircase','nSessionsHighPsych','nSessionsLowPsych','animal_ids');

figure('Position',get(0,'Screensize').*[0 0 0.5 0.5])

% a
subplot(2,3,1)
num_trial_types = size(prMean,2);
x = repmat(1:size(prMean,1),size(prMean,2),1)';
hold on
col = [1 0.7 0 ; 0.6 0.6 0.6 ; 1 0 0 ; 1 0.6 0.6];
for i = 1:num_trial_types
    errorbar(x(:,i),prMean(:,i),prSem(:,i),'k.-','CapSize',0,'MarkerSize',20,'Color',col(i,:),'LineWidth',1.5)
end
axis_tight
axis square
ylim([-0.1 1.1])
xlim([0 maxSessions+1])
xlabel('Sessions from 2P start')
ylabel('P(Lick)')
set(gca,'XTick',1:size(x,1),'YTick',[0:0.25:1])

% b
subplot(2,3,2)
hold on
line([min(x(:)) max(x(:))]+[-0.25 0.25],1*[1 1],'Color',[0.6 0.6 0.6])
col = [1 0.7 0 ; 1 0 0 ; 1 0.6 0.6];
for i = 1:num_trial_types-1
    errorbar(x(:,i),dpMean(:,i),dpSem(:,i),'k.-','CapSize',0,'MarkerSize',20,'Color',col(i,:),'LineWidth',1.5)
end
axis_tight
axis square
ylim([-0.1 4.1])
xlim([0 maxSessions+1])
xlabel('Sessions from 2P start')
ylabel('d-prime')
set(gca,'XTick',1:size(x,1))

% c
subplot(2,3,3)
col = [1 0.7 0 ; 0.6 0.6 0.6 ; 1 0 0 ; 1 0.6 0.6];
allY = [first200 last200 first100 last100];
clear dpVsZeroStats
for i = 1:size(allY,2)
    dpVsZeroStats(i) = qstat(allY(:,i)-1);
end
dpVsZeroP = [dpVsZeroStats(:).p];
idcs2Compare = [1 2 ; 3 4];
clear dpLearningStats
sameFlags = {firstLastSame200 firstLastSame100};
for i = 1:2
    tmp = allY(:,idcs2Compare(i,:));
    tmp = tmp(~sameFlags{i},:);
    dpLearningStats(i) = qstat(tmp);
end
dpLearningP = [dpLearningStats(:).p];
first100 = first100(~isnan(first100));
hold on
line([0 9],[1 1],'Color',[0.4 0.4 0.4])
scatter_errorbar(repmat([1 3],numel(first200),1),[first200 last200],'Offset',[-0.5 0.5],'MarkerSize',6,'ErrorColor',col(3,:),'LineStyle','-','LineColor',[0.6 0.6 0.6],'ScatterColor',[0.6 0.6 0.6],'ScatterEdgeColor',[0.6 0.6 0.6])
scatter_errorbar(repmat([6 8],numel(first100(~firstLastSame100)),1),[first100(~firstLastSame100) last100(~firstLastSame100)],'Offset',[-0.5 0.5],'MarkerSize',6,'ErrorColor',col(4,:),'LineStyle','-','LineColor',[0.6 0.6 0.6],'ScatterColor',[0.6 0.6 0.6],'ScatterEdgeColor',[0.6 0.6 0.6])
scatter_errorbar(repmat([6 8],numel(first100(firstLastSame100)),1),[first100(firstLastSame100) last100(firstLastSame100)],'errorbar',false,'MarkerSize',6,'LineStyle','none','LineColor',[0.6 0.6 0.6],'ScatterColor',[0.6 0.6 0.6],'ScatterEdgeColor',[0.6 0.6 0.6])
xlim([0 9])
axis_tight;
axis square
box off
set(gca,'XTick',[1 3 6 8],'XTickLabels',{'First' 'Last' 'First' 'Last'})
xtickangle(30)
plot_sig([1 3 6 8],dpVsZeroP)
plot_sig([2 7],dpLearningP)
ylabel('d-prime')
xlabel('200 cells     100 cells')

% d
subplot(2,3,4)
hold on
col = [1 0.7 0 ; 0.6 0.6 0.6 ; 1 0 0 ; 1 0.6 0.6];
for i = 1:num_trial_types
    errorbar(x(:,i),rxnAvgMean(:,i),rxnAvgSem(:,i),'k.-','CapSize',0,'MarkerSize',20,'Color',col(i,:),'LineWidth',1.5)
end
axis_tight;
ylim([-0.1 1.1])
axis square
xlim([0 maxSessions+1])
xlabel('Sessions from 2P start')
ylabel('Reaction time mean (s)')
set(gca,'XTick',1:size(x,1),'YTick',[0:0.25:1])

% e
subplot(2,3,5)
hold on
col = [1 0.7 0 ; 0.6 0.6 0.6 ; 1 0 0 ; 1 0.6 0.6];
for i = 1:num_trial_types
    errorbar(x(:,i),rxnSdMean(:,i),rxnSdSem(:,i),'k.-','CapSize',0,'MarkerSize',20,'Color',col(i,:),'LineWidth',1.5)
end
[~,yl] = axis_tight;
ylim([-0.1 1.1])
axis square
xlim([0 maxSessions+1])
xlabel('Sessions from 2P start')
ylabel('Reaction time s.d. (s)')
set(gca,'XTick',1:size(x,1),'YTick',[0:0.25:1])

% f
subplot(2,3,6)
nSessions2P = cellfun(@(x) size(x,1),all_pr);
orderIdx = nan(numel(animal_ids),1);
for a = 1:numel(animal_ids)
    idx = find(strcmp(animal_ids{a},learning_1p.animal_ids));
    if ~isempty(idx)
        orderIdx(a) = idx;
    end
end
nSessions2P = nSessions2P(~isnan(orderIdx));
orderIdx = orderIdx(~isnan(orderIdx));

y = [learning_1p.nSessionsStaircase(orderIdx) learning_1p.nSessionsHighPsych(orderIdx) learning_1p.nSessionsLowPsych(orderIdx) nSessions2P];
boxplot([y sum(y,2)],'Positions',[1 2 3 4 5],'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
boxplotEdit('Median','Width',0.5,'Median','Color',[0.7 0.7 0.7],'Median','LineWidth',2,...
    'Box','Color',[0 0 0],'Box','LineWidth',10,'Whiskers','LineWidth',2,'Whisker','Color',[0 0 0])
axis_tight
axis square
box off
set(gca,'XTickLabels',{'1P Staircase' '1P High Power psych.' '1P Low power psych.' '2P 200/100 cells' 'All'})
xtickangle(30)
ylabel('No. sessions')
