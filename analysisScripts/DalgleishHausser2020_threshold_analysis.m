%% Find thresholds
fprintf(['- Finding thresholds...\n'])

if any([~exist('tw_response_zs') ~exist('tw_response') ~exist('num_cells')])
    fprintf(['- Importing saved data...\n'])
    load('DalgleishHausser2020_imaging_raw.mat','tw_response_zs','tw_response','num_cells')
else
    fprintf(['- Using data in workspace...\n'])
end

% params
tholds = [1:0.1:3];                                                     % thresholds to sweep through
nP = 10000;                                                             % no. permutations

% setup
ng_var = num_cells==0;
pop_size = cellfun(@(x) size(x,1),tw_response_zs(:,1));                 % population size (i.e. no. ROIs)
Cr = cellfun(@(x) find(x==0),tw_response(:,ng_var),'UniformOutput',0);  % Correct reject indices for Catch trials
nCr = cellfun(@(x) sum(x==0),tw_response(:,ng_var));                    % no. correct rejects
num_sessions = size(tw_response,1);

% run permutations
nTholds = numel(tholds);
try
    load('DalgleishHausser2020_permutedFPRate.mat')
catch
    [e_fp_perm,i_fp_perm] = deal(zeros(num_sessions,nP,nTholds));
    rng('default')
    for t = 1:nTholds
        for s = 1:num_sessions
            for p = 1:nP
                % train/test split
                [trainInd,valInd] = dividerand(nCr(s),0.8,0.2);
                trainInd = Cr{s}(trainInd);
                valInd = Cr{s}(valInd);
                
                % "train" cell thresholds using this scale factor
                mu = mean(tw_response_zs{s,ng_var}(:,trainInd),2);
                sd = std(tw_response_zs{s,ng_var}(:,trainInd),[],2);
                pos_thold_perm = mu + tholds(t) * sd;
                neg_thold_perm = mu - tholds(t) * sd;
                
                % "test" cell thresholds using this scale factor
                e_fp_perm(s,p,t) = mean(sum(tw_response_zs{s,ng_var}(:,valInd) > pos_thold_perm,1),2);
                i_fp_perm(s,p,t) = mean(sum(tw_response_zs{s,ng_var}(:,valInd) < neg_thold_perm,1),2);
            end
        end
    end
    % convert to proportion of cells
    e_fp_perm = e_fp_perm./pop_size;
    i_fp_perm = i_fp_perm./pop_size;
    %save([backdir(which('DalgleishHausser2020_imaging_raw.mat'),1) filesep 'DalgleishHausser2020_permutedFPRate.mat'],'e_fp_perm','i_fp_perm','tholds','nP')
end

%% Define and plot cross-validated thresholds (scaling factor with median 5% FP across shuffles)
if ~exist('num_sessions')
    num_sessions = size(e_fp_perm,1);
end
desiredFpRate = 0.05; % FP rate you want to aim for (this is used to find the x-value corresponding to y = 0.05 point of fitted curves below)

% use median fp rate across each permutation
eFP = squeeze(median(e_fp_perm,2)); 
iFP = squeeze(median(i_fp_perm,2));

% fit curves
xFit = linspace(min(tholds(:)),max(tholds(:)),1000);
[eThold,iThold] = deal(zeros(num_sessions,1));
[eCurves,iCurves] = deal(zeros(num_sessions,numel(xFit)));

for s = 1:num_sessions    
    % fitted curves
    eCurves(s,:) = interp1(tholds(:),eFP(s,:)',xFit,'pchip');
    iCurves(s,:) = interp1(tholds(:),iFP(s,:)',xFit,'pchip');
    
    % find scaling factors
    eThold(s) = fzero(@(xi) interp1(tholds(:),eFP(s,:)',xi,'pchip')-desiredFpRate,2);
    iThold(s) = fzero(@(xi) interp1(tholds(:),iFP(s,:)',xi,'pchip')-desiredFpRate,2);
end

% return actual thresholds for each neuron
[pos_thold,neg_thold] = deal(cell(num_sessions,1));
for s = 1:num_sessions
    cr = tw_response{s,ng_var}==0;
    mu = mean(tw_response_zs{s,ng_var}(:,cr),2);
    sd = std(tw_response_zs{s,ng_var}(:,cr),[],2);
    pos_thold{s} = mu + eThold(s) * sd;
    neg_thold{s} = mu - iThold(s) * sd;
end

% plot
figure
subplot(1,2,1)
hold on
xRaw = repmat(tholds,num_sessions,1);
x = repmat(xFit,size(eCurves,1),1);
plot(x',eCurves','Color',[1 0.6 0.6])
scatter(xRaw(:),eFP(:),150,'r.')
[xl,yl] = axis_tight;
yt = get(gca,'YTick');
boxplot(eThold,'orientation','horizontal','Positions',min(eFP(:)),'Widths',0.005,'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
xlim(xl); ylim(yl)
axis square
line(xl,desiredFpRate * [1 1],'Color',[0.4 0.4 0.4])
scatter(eThold,desiredFpRate*(1+0*eThold),190,'k.')
box off
set(gca,'XTick',[1:0.5:3],'YTick',yt(1):0.05:yt(end),'YTickLabels',yt(1):0.05:yt(end))
xlabel('Threshold S.D. scale factor')
ylabel('Median FP rate')
title('Activation threshold')

subplot(1,2,2)
hold on
xRaw = repmat(tholds,num_sessions,1);
x = repmat(xFit,size(iCurves,1),1);
plot(x',iCurves','Color',[0.6 0.6 1])
scatter(xRaw(:),iFP(:),150,'b.')
[xl,yl] = axis_tight;
ylim([0 0.2])
yt = get(gca,'YTick');
boxplot(iThold,'orientation','horizontal','Positions',min(iFP(:)),'Widths',0.005,'BoxStyle','filled','MedianStyle','line','Colors',[0.4 0.4 0.4])
xlim(xl); ylim(yl)
axis square
line(xl,desiredFpRate * [1 1],'Color',[0.4 0.4 0.4])
scatter(iThold,desiredFpRate*(1+0*iThold),190,'k.')
box off
set(gca,'XTick',[1:0.5:3],'YTick',yt(1):0.05:yt(end),'YTickLabels',yt(1):0.05:yt(end))
xlabel('Threshold S.D. scale factor')
ylabel('Median FP rate')
title('Suppression threshold')

%% Define and plot non-cross-validated thresholds
[tmp_e_fp,tmp_i_fp] = deal(zeros(num_sessions,nTholds));
for t = 1:nTholds
    for s = 1:num_sessions
        mu = mean(tw_response_zs{s,ng_var}(:,Cr{s}),2);
        sd = std(tw_response_zs{s,ng_var}(:,Cr{s}),[],2);
        pos_thold_perm = mu + tholds(t) * sd;
        neg_thold_perm = mu - tholds(t) * sd;
        
        tmp_e_fp(s,t) = mean(sum(tw_response_zs{s,ng_var}(:,Cr{s}) > pos_thold_perm,1),2);
        tmp_i_fp(s,t) = mean(sum(tw_response_zs{s,ng_var}(:,Cr{s}) < neg_thold_perm,1),2);
    end
end
tmp_e_fp = tmp_e_fp./pop_size;
tmp_i_fp = tmp_i_fp./pop_size;

[eThold_noCV,iThold_noCV] = deal(zeros(num_sessions,1));
xFit = linspace(min(tholds(:)),max(tholds(:)),1000);
[eCurves_noCV,iCurves_noCV] = deal(zeros(num_sessions,numel(xFit)));
for s = 1:num_sessions    
    % find thresholds
    eCurves_noCV(s,:) = interp1(tholds(:),tmp_e_fp(s,:)',xFit,'pchip');
    iCurves_noCV(s,:) = interp1(tholds(:),tmp_i_fp(s,:)',xFit,'pchip');
    
    eThold_noCV(s) = fzero(@(xi) interp1(tholds(:),tmp_e_fp(s,:)',xi,'pchip')-0.05,2);
    iThold_noCV(s) = fzero(@(xi) interp1(tholds(:),tmp_i_fp(s,:)',xi,'pchip')-0.05,2);
end

figure

subplot(1,2,1)
hold on
x = repmat(xFit,num_sessions,1);
plot(x',eCurves_noCV','Color',[1 0 0])
axis square
[xl,~] = axis_tight;
line(xl,0.05*[1 1],'Color',[0.6 0.6 0.6])
box off
xlabel('Threshold S.D. scale factor')
ylabel('Proportion responsive neurons')
title('Activation threshold')

subplot(1,2,2)
hold on
plot(x',iCurves_noCV','Color',[0 0 1])
axis square
[xl,~] = axis_tight;
line(xl,0.05*[1 1],'Color',[0.6 0.6 0.6])
box off
xlabel('Threshold S.D. scale factor')
ylabel('Proportion responsive neurons')
title('Activation threshold')

%% compare cross-validated/non-cross-validated thresholds
figure
hold on
scatter(eThold_noCV(:),eThold,300,'r.')
scatter(iThold_noCV(:),iThold,300,'b.')
l = equal_axes;
axis square
line(l,l,'Color',[0.6 0.6 0.6],'LineStyle','--')
xlabel('S.D. scaling factor (not cross-validated)')
ylabel('S.D. scaling factor (cross-validated)')
