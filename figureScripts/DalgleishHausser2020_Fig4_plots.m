%% Figure 4
if ~exist('isRun') || (exist('isRun') && ~isfield(isRun,'figure2'))
    run('DalgleishHausser2020_Fig2_3_4_analysis')
end

plotOptions             = struct;
plotOptions.plotThresh  = false;
plotOptions.CIthresh    = false; 
plotOptions.dataColor   = [0.6 0.6 0.6];

scrsz = get(0,'Screensize');
figure('Position',[scrsz(3)*0.12 0 scrsz(3:end).*[0.6 0.3]])

% a,b,c
xlabels = {'Target activation' 'Background activation' 'Background suppression'};
for i = 1:nPredictors
    subplot(1,nPredictors+1,i)
    plotPsych(predictorFits(i),plotOptions)
    title(['R2 = ' num2str(round(r2(i),2))])
    axis_tight
    axis square
    xlabel(xlabels{i})
    ylabel('P(lick)')
    set(gca,'YTick',[0:0.25:1])
end

 % d
subplot(1,nPredictors+1,i+1)
r2SigVal = permTestP(indModel_r2PermVal,0);
r2SigFit = permTestP(indModel_r2PermFit,0);
hold on
line([0 5],[0 0],'Color',[0.6 0.6 0.6])
xVals = [1 2.5 4];
offset = 0.3;
bp = boxplot(indModel_r2PermFit,'BoxStyle','filled','Widths',0.5,'MedianStyle','line','Colors',[0.4 0.4 0.4],'Symbol','','Positions',xVals-offset);
boxplotEdit('Box','Color',[0.6 0.6 0.6],'Box','LineWidth',8,'Median','LineWidth',3,...
    'Median','Color',[0.6 0.6 0.6],'Median','Width',0.45,...
    'Whiskers','LineWidth',2,'Whikers','Color',[0.6 0.6 0.6]);
bp = boxplot(indModel_r2PermVal,'BoxStyle','filled','Widths',0.5,'MedianStyle','line','Colors',[0.4 0.4 0.4],'Symbol','','Positions',xVals+offset);
boxplotEdit('Box','Color',[0 0 0],'Box','LineWidth',8,'Median','LineWidth',3,...
    'Median','Color',[0 0 0],'Median','Width',0.45,...
    'Whiskers','LineWidth',2,'Whikers','Color',[0.6 0.6 0.6]);
xlim([0 5])
set(gca,'XTick',xVals,'XTickLabels',{'Target activation' 'Network activation' 'Network suppression'})
xtickangle(30)
axis square
box off
ylabel('Test R2')
ylim([-1.1 1.1])
plot_sig(xVals-offset,r2SigFit)
plot_sig(xVals+offset,r2SigVal)
