%% Figure 3
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

scrsz = get(0,'Screensize');
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.6 0.5]])

sz = 20;
exclude_few = true;
corr_col = [0 0.5 0.8 ; 0.8 0.1 0.1];

x = log(e_t_n(:,~ng_var));
y_e_all = e_all_matched(:,~ng_var);
y_i_all = i_all_matched(:,~ng_var);
y_t_e = e_t_matched(:,~ng_var);
y_i = i_off_matched(:,~ng_var);
y_e = e_off_matched(:,~ng_var);
y_d = e_all_matched ./ i_all_matched;
y_d_catch = y_d(:,ng_var);
y_d = y_d(:,~ng_var);
flag = isnan(y_i) | isnan(y_e) | isnan(y_d) | isnan(x);
if exclude_few
    flag = flag | e_t_n(:,~ng_var)<1;
    x(flag) = nan; y_e(flag) = nan; y_e(flag) = nan; y_d(flag) = nan;
end

% a
subplot(2,3,1)
hold on
shadedErrorBar([min(x(:)) max(x(:))],nanmean(e_all_matched(:,ng_var))*[1 1],sem(e_all_matched(:,ng_var))*[1 1],{'Color',[1 0.8 0.8]})
scatter(x(:),y_e_all(:),100,'k.','MarkerEdgeColor',[0.6 0.6 0.6])
errorbar(nanmean(x,1),nanmean(y_e_all,1),sem(y_e_all,1),sem(y_e_all,1),sem(x,1),sem(x,1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(2,:))
eFitAll = lin_fit_plot(x(~flag),y_e_all(~flag),'Scatter',0);
[xl,~] = axis_tight;
axis square
set(gca,'XTick',log([0.1 1 10 100]),'XTickLabels',[0.1 1 10 100])
xlabel('No. activated target neurons (log)')
ylabel({'P(Activated) Lick matched' 'All neurons'})
yTmp = [nanmean(y_e_all,2) e_all_matched(:,ng_var)];
yTmp(any(isnan(yTmp),2),:) = [];
errorbar(max(xl),nanmean(yTmp(:,2),1),sem(yTmp(:,2),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[1 0.8 0.8])
errorbar(max(xl),nanmean(yTmp(:,1),1),sem(yTmp(:,1),1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(2,:))
title_stats(qstat(yTmp))

% b
subplot(2,3,2)
hold on
shadedErrorBar([min(x(:)) max(x(:))],nanmean(i_all_matched(:,ng_var))*[1 1],sem(i_all_matched(:,ng_var))*[1 1],{'Color',[0.8 0.8 1]})
scatter(x(:),y_i_all(:),100,'k.','MarkerEdgeColor',[0.6 0.6 0.6])
errorbar(nanmean(x,1),nanmean(y_i_all,1),sem(y_i_all,1),sem(y_i_all,1),sem(x,1),sem(x,1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(1,:))
iFitAll = lin_fit_plot(x(~flag),y_i_all(~flag),'Scatter',0);
axis_tight;
axis square
set(gca,'XTick',log([0.1 1 10 100]),'XTickLabels',[0.1 1 10 100])
xlabel('No. activated target neurons (log)')
ylabel({'P(Suppressed) Lick matched' 'All neurons'})
yTmp = [nanmean(y_i_all,2) i_all_matched(:,ng_var)];
yTmp(any(isnan(yTmp),2),:) = [];
errorbar(max(xl),nanmean(yTmp(:,2),1),sem(yTmp(:,2),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[0.8 0.8 1])
errorbar(max(xl),nanmean(yTmp(:,1),1),sem(yTmp(:,1),1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(1,:))
title_stats(qstat(yTmp))

% c
subplot(2,3,3)
hold on
shadedErrorBar([min(x(:)) max(x(:))],nanmean(y_d_catch)*[1 1],sem(y_d_catch)*[1 1],{'Color',[0.7 0.7 0.7]})
scatter(x(:),y_d(:),100,'k.','MarkerEdgeColor',[0.6 0.6 0.6])
errorbar(nanmean(x,1),nanmean(y_d,1),sem(y_d,1),sem(y_d,1),sem(x,1),sem(x,1),'k.','MarkerSize',sz,'CapSize',0,'Color',[0 0 0])
dFit = lin_fit_plot(x(~flag),y_d(~flag),'Scatter',0);
[xl,~] = axis_tight;
axis square
set(gca,'XTick',log([0.1 1 10 100]),'XTickLabels',[0.1 1 10 100])
xlabel('No. activated target neurons (log)')
ylabel({'Activation / Suppression' 'All neurons'})
yTmp = [nanmean(y_d,2) y_d_catch];
yTmp(any(isnan(yTmp),2),:) = [];
errorbar(max(xl),nanmean(yTmp(:,2),1),sem(yTmp(:,2),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[0.7 0.7 0.7])
errorbar(max(xl),nanmean(yTmp(:,1),1),sem(yTmp(:,1),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[0 0 0])
title_stats(qstat(yTmp))

% d
subplot(2,3,4)
hold on
shadedErrorBar([min(x(:)) max(x(:))],nanmean(e_off_matched(:,ng_var))*[1 1],sem(e_off_matched(:,ng_var))*[1 1],{'Color',[1 0.8 0.8]})
scatter(x(:),y_e(:),100,'k.','MarkerEdgeColor',[0.6 0.6 0.6])
errorbar(nanmean(x,1),nanmean(y_e,1),sem(y_e,1),sem(y_e,1),sem(x,1),sem(x,1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(2,:))
eFit = lin_fit_plot(x(~flag),y_e(~flag),'Scatter',0);
[xl,~] = axis_tight;
axis square
set(gca,'XTick',log([0.1 1 10 100]),'XTickLabels',[0.1 1 10 100])
xlabel('No. activated target neurons (log)')
ylabel({'P(Activated) Lick matched' 'Network'})
yTmp = [nanmean(y_e,2) e_off_matched(:,ng_var)];
yTmp(any(isnan(yTmp),2),:) = [];
errorbar(max(xl),nanmean(yTmp(:,2),1),sem(yTmp(:,2),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[1 0.8 0.8])
errorbar(max(xl),nanmean(yTmp(:,1),1),sem(yTmp(:,1),1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(2,:))
title_stats(qstat(yTmp))

% e
subplot(2,3,5)
hold on
shadedErrorBar([min(x(:)) max(x(:))],nanmean(i_off_matched(:,ng_var))*[1 1],sem(i_off_matched(:,ng_var))*[1 1],{'Color',[0.8 0.8 1]})
scatter(x(:),y_i(:),100,'k.','MarkerEdgeColor',[0.6 0.6 0.6])
errorbar(nanmean(x,1),nanmean(y_i,1),sem(y_i,1),sem(y_i,1),sem(x,1),sem(x,1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(1,:))
iFit = lin_fit_plot(x(~flag),y_i(~flag),'Scatter',0);
axis_tight;
axis square
set(gca,'XTick',log([0.1 1 10 100]),'XTickLabels',[0.1 1 10 100])
xlabel('No. activated target neurons (log)')
ylabel({'P(Suppressed) Lick matched' 'Network'})
yTmp = [nanmean(y_i,2) i_off_matched(:,ng_var)];
yTmp(any(isnan(yTmp),2),:) = [];
errorbar(max(xl),nanmean(yTmp(:,2),1),sem(yTmp(:,2),1),'k.','MarkerSize',sz,'CapSize',0,'Color',[0.8 0.8 1])
errorbar(max(xl),nanmean(yTmp(:,1),1),sem(yTmp(:,1),1),'k.','MarkerSize',sz,'CapSize',0,'Color',corr_col(1,:))
title_stats(qstat(yTmp))

%% Figure 3 - figure supplement 1a-h
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

c = [0.2 0.7 0.1 ; 0.6 0.6 0.6];

scrsz = get(0,'Screensize');
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.6 0.6]])

% a
h = subplot(2,4,1);
hold on
x = e_t_n(:,~ng_var);
xMean = nanmean(x,1);
xErr = sem(x,1);
xCatch = [xMean(end)-xErr(end) xMean(1)+xErr(1)];
shadedErrorBar(xCatch,nanmean(e_off_hit(:,ng_var))*[1 1],sem(e_off_hit(:,ng_var))*[1 1],{'Color',c(1,:)});
shadedErrorBar(xCatch,nanmean(e_off_miss(:,ng_var))*[1 1],sem(e_off_miss(:,ng_var))*[1 1],{'Color',c(2,:)});
yErr = sem(e_off_hit(:,~ng_var),[],1);
errorbar(xMean,nanmean(e_off_hit(:,~ng_var),1),yErr,yErr,xErr,xErr,'ko-','Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',6,'CapSize',0,'LineWidth',1.5);
yErr = sem(e_off_miss(:,~ng_var),[],1);
errorbar(xMean,nanmean(e_off_miss(:,~ng_var),1),yErr,yErr,xErr,xErr,'ko-','Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',6,'CapSize',0,'LineWidth',1.5);
axis square
box off
xlabel({'No. activated' 'target neurons'})
ylabel({'P(activated)' 'background'})
axis_tight
set(gca,'XScale','log','XTick',[1 10 100],'XTickLabels',[1 10 100])
axis_tight

% b
subplot(2,4,2)
hold on
y = [e_off_hit(:,thold_var) e_off_miss(:,thold_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c,'ScatterColor',c,'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'50 neuron hit' '50 cell miss'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(activated)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

% c
subplot(2,4,3)
hold on
y = [e_off_hit(:,ng_var) e_off_miss(:,ng_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c,'ScatterColor',c,'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'Catch FA' 'Catch CR'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(activated)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

% d
subplot(2,4,4)
hold on
y = [e_off_hit(:,thold_var) e_off_hit(:,ng_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c(1,:),'ScatterColor',c(1,:),'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'50 neuron hit' 'Catch FA'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(activated)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

% e
subplot(2,4,5)
hold on
x = e_t_n(:,~ng_var);
xMean = nanmean(x,1);
xErr = sem(x,1);
xCatch = [xMean(end)-xErr(end) xMean(1)+xErr(1)];
shadedErrorBar(xCatch,nanmean(i_off_hit(:,ng_var))*[1 1],sem(i_off_hit(:,ng_var))*[1 1],{'Color',c(1,:)});
shadedErrorBar(xCatch,nanmean(i_off_miss(:,ng_var))*[1 1],sem(i_off_miss(:,ng_var))*[1 1],{'Color',c(2,:)});
yErr = sem(e_off_hit(:,~ng_var),[],1);
errorbar(xMean,nanmean(i_off_hit(:,~ng_var),1),yErr,yErr,xErr,xErr,'ko-','Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',6,'CapSize',0,'LineWidth',1.5);
yErr = sem(e_off_miss(:,~ng_var),[],1);
errorbar(xMean,nanmean(i_off_miss(:,~ng_var),1),yErr,yErr,xErr,xErr,'ko-','Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',6,'CapSize',0,'LineWidth',1.5);
axis square
box off
xlabel({'No. activated' 'target neurons'})
ylabel({'P(suppressed)' 'background'})
axis_tight
set(gca,'XScale','log','XTick',[1 10 100],'XTickLabels',[1 10 100])
axis_tight

% f
subplot(2,4,6)
hold on
y = [i_off_hit(:,thold_var) i_off_miss(:,thold_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c,'ScatterColor',c,'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'50 neuron hit' '50 neuron miss'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(suppressed)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

% g
h = subplot(2,4,7);
hold on
y = [i_off_hit(:,ng_var) i_off_miss(:,ng_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c,'ScatterColor',c,'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'Catch FA' 'Catch CR'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(suppressed)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

% h
subplot(2,4,8)
hold on
y = [i_off_hit(:,thold_var) i_off_hit(:,ng_var)];
y(any(isnan(y),2),:) = [];
scatter_errorbar([],y,'LineWidth',1,'LineColor',[0.6 0.6 0.6],'ErrorColor',c(1,:),'ScatterColor',c(1,:),'Offset',[-0.2 0.2],'MarkerSize',7);
set(gca,'XTick',[1 2],'XTickLabels',{'50 neuron hit' 'Catch FA'})
title_stats(qstat(y));
xlim([0.5 2.5])
xtickangle(30)
ylabel({'P(suppressed)' 'background'})
ax = get(gca);
ax.YAxis.Exponent = -2;
axis square
plot_sig([1 2],qstat(y))

%% Figure 3 - figure supplement 1i - n
% NB since this plot uses "raw traces", please download raw data from:
%
% https://doi.org/10.6084/m9.figshare.13128950
%
% and import using DalgleishHausser2020_importProcessing.m in this repo

if exist('import_raw') && import_raw
    try
        if ~exist('spont_lick_cc')
            run('DalgleishHausser2020_importProcessing.m')
        end
        
        scrsz = get(0,'Screensize');
        figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.7 0.45]])
        subplot(2,4,[1 2])
        epoch = [4700:5700];
        hold on
        s = 5;
        idx = 1;
        [vals,order] = sort(spont_lick_cc{s},'Descend');
        [~,f] = bl_normalize_trace(session(s).f_sub(order(idx),:),running_baseline_s,imaging_rate_rounded);
        f = zscore(f,[],2);
        f = gaussfilt1d(f(:,~stim_trace{s}),imaging_rate_rounded*1);
        l = session(s).lick_trace(~stim_trace{s});
        if ischar(epoch)
            plot(normalise(f))
            scatter(find(l),0*find(l),'ko')
        else
            f = f(epoch);
            plot(f)
            lt = find(l(epoch));
            line([lt ; lt],0*[lt ; lt] + min(f) + [-0.3 ; 0],'Color',[1 0 0])
        end
        axis_tight
        box off
        set(gca,'XTick',[0 imaging_rate*10],'XTickLabels',[0 10],'YTick',[0 2])
        ylabel('dF/sF')
        
        
        subplot(2,4,5)
        y = cell2mat(spont_lick_cc);
        lims = dat2lims(y);
        histogram(y,min(y):0.01:max(y),'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',1)
        axis square
        axis_tight
        box off
        xlabel('Correlation with spont. licks')
        set(gca,'YScale','log')
        ylabel('Count (log)')
        
        subplot(2,4,6)
        hold on
        y = cell2mat(spont_lick_p);
        histogram(y,[0:0.05:1],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',1)
        axis square
        axis_tight
        box off
        %set(gca,'XScale','log','YScale','log')
        line(0.05*[1 1],get(gca,'YLim'),'Color',[1 0 0])
        ylim([0 max(get(gca,'YLim'))])
        xlabel('Spont. lick correlation p-value')
        ylabel('Count')
        
        numBins = 10;
        CCbinnedTraces = cell(num_sessions,numBins);
        binnedCC = [];
        for s = 1:num_sessions
            [binnedCC(s,:),~,~,~,yVals] = prctileYwithX(spont_lick_cc{s},1:numel(spont_lick_cc{s}),numBins);
            for b = 1:numBins
                CCbinnedTraces{s,b} = mean_lick_traces_zs{s}(yVals{b},:);
            end
        end
        binnedCCMean = mean(binnedCC,1);
        
        subplot(2,4,3)
        cMap = parula(numBins + 2);
        cMap = cMap(2:end-1,:);
        for i = 1:numBins
            hold on
            y = cell2mat(CCbinnedTraces(:,i));
            shadedErrorBar(1:size(y,2),mean(y,1),sem(y,1),{'Color',cMap(i,:)})
        end
        line(pre_frames_sta.ps+1 * [1 1],get(gca,'YLim'),'Color',[1 0 0])
        axis_tight
        axis square
        set(gca,'XTick',[0 imaging_rate],'XTickLabels',[0 1])
        cb = colorbar;
        tl = round(min(binnedCCMean),2):0.01:round(max(binnedCCMean),2);
        cb.Ticks = linspace(min(cb.Ticks),max(cb.Ticks),numel(tl));
        cb.TickLabels = tl;
        cb.Location = 'east';
        
        subplot(2,4,7)
        y = cell2mat(tw_spont_lick_traces);
        shadedErrorBar(1:size(y,2),nanmean(y,1),sem(y,1),{'Color',[0 0 0]})
        axis_tight
        axis square
        line(pre_frames_sta.ps+1 * [1 1],get(gca,'YLim'),'Color',[1 0 0])
        box off
        
        [crossCorrDffPosSig,crossCorrDffNegSig,crossCorrLagPosSig,crossCorrLagNegSig,...
            spont_lick_cc_sig_pos,spont_lick_cc_sig_neg] = deal({});
        mid = ceil(size(crossCorrDff{1},2)/2);
        epoch = round(mid + [-imaging_rate*2:imaging_rate*2]);
        for s = 1:num_sessions
            crossCorrDffPosSig{s,1} = crossCorrDff{s}(spont_lick_cc{s}>0 & lick_exclude{s},:);
            crossCorrDffNegSig{s,1} = crossCorrDff{s}(spont_lick_cc{s}<0 & lick_exclude{s},:);
            
            spont_lick_cc_sig_pos{s,1} = spont_lick_cc{s}(spont_lick_cc{s}>0 & lick_exclude{s});
            spont_lick_cc_sig_neg{s,1} = spont_lick_cc{s}(spont_lick_cc{s}<0 & lick_exclude{s});
            
            [~,crossCorrLagPosSig{s,1}] = max(crossCorrDffPosSig{s,1}(:,epoch),[],2);
            crossCorrLagPosSig{s,1} = (crossCorrLagPosSig{s,1} + min(epoch) - mid)/imaging_rate;
            [~,crossCorrLagNegSig{s,1}] = min(crossCorrDffNegSig{s,1}(:,epoch),[],2);
            crossCorrLagNegSig{s,1} = (crossCorrLagNegSig{s,1} + min(epoch) - mid)/imaging_rate;
        end
        
        subplot(2,4,4)
        hold on
        y = cell2mat(crossCorrDffPosSig);
        x = repmat(lags/imaging_rate,size(y,1),1);
        shadedErrorBar(mean(x,1),mean(y,1),sem(y,1),{'Color',cMap(end,:)})
        y = cell2mat(crossCorrDffNegSig);
        x = repmat(lags/imaging_rate,size(y,1),1);
        shadedErrorBar(mean(x,1),mean(y,1),sem(y,1),{'Color',cMap(1,:)})
        axis_tight
        axis square
        box off
        line([0 0],get(gca,'YLim'))
        xlabel('Lag (s)')
        ylabel({'Lick - dF/F' 'cross-correlation (norm.)'})
        
        subplot(2,4,8)
        hold on
        yPos = cell2mat(crossCorrLagPosSig);
        yNeg = cell2mat(crossCorrLagNegSig);
        histogram(yPos,200,'Normalization','cdf','FaceColor','none','DisplayStyle','stairs','EdgeColor',cMap(end,:),'LineWidth',1.5)
        histogram(yNeg,200,'Normalization','cdf','FaceColor','none','DisplayStyle','stairs','EdgeColor',cMap(1,:),'LineWidth',1.5)
        axis_tight
        axis square
        box off
        title_stats(qstat({yPos yNeg}))
        xlabel('XCorr peak time (s)')
        ylabel('Percentile')
    catch
        fprintf(['\nThese panels require raw data, please download from:\n'...
            '\nhttps://doi.org/10.6084/m9.figshare.13128950\n'...
            '\nIf you have downloaded this data, ensure the correct data directory is referenced in DalgleishHausser2020_importProcessing.m'...
            '\nand set "import_raw" to true\n']);
    end
else
    fprintf(['\nThese panels require raw data, please download from:\n'...
            '\nhttps://doi.org/10.6084/m9.figshare.13128950\n'...
            '\nIf you have downloaded this data, ensure the correct data directory is referenced in DalgleishHausser2020_importProcessing.m'...
            '\nand set "import_raw" to true\n']);
end

%% Figure 3 - figure supplement 3
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

% xy vs p(response) plot
posLims     = [-0.1 0.1];
negLims     = [-0.05 0.05];
posCmap     = makeColormap([[0 0 0];[1 1 1];[1 0.3 0.3]],999);
negCmap     = makeColormap([[0 0 0];[1 1 1];[0.3 0.3 1]],999);

scrsz = get(0,'Screensize');
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.6 0.5]])

% E heatmap
ax = [];
ax(1) = subplot(2,3,1);
im = nanmean(xyzE,3);
imagesc(xyBins,zDepths,im,posLims);
yt = get(gca,'YTick');
colormap(ax(1),posCmap)
axis square
xl = get(gca,'XLim');
yl = get(gca,'YLim');
wid = 100-max(yl);
hold on
line([5 5],yl,'Color',[0.6 0.6 0.6])
xlim(xl)
ylim(yl)
xlabel({'Lateral distance from' 'nearest light (µm)'})
ylabel({'Axial distance from' 'nearest light (µm)'})
cb = colorbar('Location','southoutside','Orientation','horizontal','Ticks',round(sort([posLims 0]),2));
cb.Position = cb.Position .* [1 0.75 0.5 0.75];
cb.Label.String = {'P(activated)' 'spont sub.'};
title('Activation')
set(gca,'XTick',[0:50:100]-5,'XTickLabels',[0:50:100],'YTick',[0:50:100],'YDir','normal')

% I heatmap
ax(2) = subplot(2,3,2);
im = nanmean(xyzI,3);
imagesc(xyBins,zDepths,im,negLims)
colormap(ax(2),negCmap)
axis square
xl = get(gca,'XLim');
yl = get(gca,'YLim');
wid = 100-max(yl);
hold on
line([5 5],yl,'Color',[0.6 0.6 0.6])
xlim(xl)
ylim(yl)
xlabel({'Lateral distance from' 'nearest light (µm)'})
ylabel({'Axial distance from' 'nearest light (µm)'})
cb = colorbar('Location','southoutside','Orientation','horizontal','Ticks',round(sort([negLims 0]),2));
cb.Position = cb.Position .* [1 0.75 0.5 0.75];
cb.Label.String = {'P(suppressed)' 'spont sub.'};
title('Suppression')
set(gca,'XTick',[0:50:100]-5,'XTickLabels',[0:50:100],'YTick',[0:50:100],'YDir','normal')

% lateral distance vs activation
subplot(2,3,3)
yLim = [-0.02 0.23];
yE = squeeze(nanmean(xyE,1))';
yI = -squeeze(nanmean(xyI,1))';
[pE,pI] = deal(0*yE(:,1));
for i = 1:size(yE,2)
    pE(i) = signrank(yE(:,i));
    pI(i) = signrank(yI(:,i));
end
pE(isnan(pE)) = []; pI(isnan(pI)) = [];
pE = pE * numel(pE);
pI = pI * numel(pI);
hold on
bar(xyBins(pE<0.05),nanmean(yE(:,pE<0.05),1),'FaceColor',[1 0.6 0.6],'EdgeColor','none')
bar(xyBins(pE>=0.05),nanmean(yE(:,pE>=0.05),1),'EdgeColor',[1 0.6 0.6],'FaceColor','none')
errorbar(xyBins,nanmean(yE,1),[],sem(yE,1),'k-','Color',[1 0.6 0.6],'LineStyle','none','CapSize',0,'LineWidth',2)
bar(xyBins(pI<0.05),nanmean(yI(:,pI<0.05),1),'FaceColor',[0.6 0.6 1],'EdgeColor','none')
bar(xyBins(pI>=0.05),nanmean(yI(:,pI>=0.05),1),'EdgeColor',[0.6 0.6 1],'FaceColor','none')
errorbar(xyBins,nanmean(yI,1),sem(yI,1),'k-','Color',[0.6 0.6 1],'LineStyle','none','CapSize',0,'LineWidth',2)
axis square
set(gca,'XTick',0:25:100)
[~,yl] = axis_tight;
ylim(yLim)
xlim([-5 105])
line([10 10]-diff(xyBins(1:2))/2,yl,'Color',[0.4 0.4 0.4])
ylabel({'Cell response probability' 'Spont. subtracted'})
xlabel({'Lateral distance from' 'nearest light (µm)'})
plot_sig(xyBins,pE)
plot_sig(xyBins,pI,0.01)
title('Across all depths')

