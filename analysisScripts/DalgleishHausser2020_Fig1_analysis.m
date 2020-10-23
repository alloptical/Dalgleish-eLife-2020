%% Figure 1 and Figure 1 - figure supplement 5
% 2P/1P behaviour

load('DalgleishHausser2020_All2PLearning_raw.mat')

% analyse example behaviour session
params.response_window     = 1.15;
params.sated_thold         = 0.7;
params.min_dur             = 0.5;
params.min_lick_t          = 0.15;
params.remove_autoreward   = false;
params.remove_ambiguous    = true;
bhv_file = 'H527_20190206_165853.mat';
[bhv_eg] = pyb_load_bhv_file(bhv_file,params);

% stims are ordered by stim_type/stim_var
% this corresponds to [1P Catch 2P_200 2P_100]
all_stim_vars = [5 6 7 7; 0 0 0 1];
num_animals = numel(all_bhv);
[all_pr,all_dp,all_rxn_avg,all_rxn_sd,pAutoReward] = deal(cell(numel(all_bhv),1));
finalNonAutoSession = nan(numel(all_bhv),1);
[lastPr,lastRxnMean,lastRxnSd] = deal(nan(num_animals,size(all_stim_vars,2)));
num_trial_types = size(all_stim_vars,2);
lastLearntIdx = nan(num_animals,1);

for a = 1:num_animals
    [all_pr{a},all_dp{a},all_rxn_avg{a},all_rxn_sd{a}] = deal(nan(numel(all_bhv(a).bhv_daily),size(all_stim_vars,2)));
    for b = 1:numel(all_bhv(a).bhv_daily)
        ng = pyb_return_trials(all_bhv(a).bhv_daily{b},6,0);
        for s = 1:size(all_stim_vars,2)
            o = pyb_return_trials(all_bhv(a).bhv_daily{b},all_stim_vars(1,s),all_stim_vars(2,s));
            if numel(o) > 10
                all_pr{a}(b,s) = mean([o(:).responded]);
                all_dp{a}(b,s) = dprime(sum([o(:).responded]),numel(o),sum([ng(:).responded]),numel(ng));
                
                if sum([o(:).responded]) >= 3
                    all_rxn_avg{a}(b,s) = nanmean([o(:).rxntime]);
                    all_rxn_sd{a}(b,s) = nanstd([o(:).rxntime]);
                end
            end
        end
    end
    
    % find last 1P/2P session where animal has learnt
    lastLearntIdx(a)           = find(all_dp{a}(:,3)>1 & ~any(isnan(all_dp{a}(:,1:3)),2),1,'Last');
    if ~isempty(lastLearntIdx(a))
        lastPr(a,:)         = all_pr{a}(lastLearntIdx(a),:);
        lastRxnMean(a,:)    = all_rxn_avg{a}(lastLearntIdx(a),:);
        lastRxnSd(a,:)      = all_rxn_sd{a}(lastLearntIdx(a),:);
    end
end

% align sessions to first one
numSessions = cellfun(@(x) size(x,1),all_pr);
maxSessions = max(numSessions);
[allPrStartAligned,allRxnAvgStartAligned,allRxnSdStartAligned] = deal(nan(maxSessions,num_trial_types,num_animals));
allDpStartAligned = nan(maxSessions,num_trial_types-1,num_animals);

for a = 1:num_animals
    thisN = size(all_pr{a},1);
    allPrStartAligned(1:thisN,:,a) = all_pr{a};
    allDpStartAligned(1:thisN,:,a) = all_dp{a}(:,[1 3:end]);
    allRxnAvgStartAligned(1:thisN,:,a) = all_rxn_avg{a};
    allRxnSdStartAligned(1:thisN,:,a) = all_rxn_sd{a};
end
% find session no. with >=3 animals (for plotting)
nAnimals    = sum(any(~isnan(allPrStartAligned),2),3);
maxSessions = find(nAnimals>=3,1,'Last');
% find time to learn 2P
tmp = squeeze(allDpStartAligned(:,2,:));
learnTime2P200 = nan(num_animals,1);
for i = 1:size(tmp,2)
    learnTime2P200(i) = find(tmp(:,i)>1,1,'First');
end
% find no. 2P sessions
nSessions2P = sum(squeeze(any(~isnan(allPrStartAligned),2)),1);

prMean      = nanmean(allPrStartAligned(nAnimals>=3,:,:),3);
prSem       = sem(allPrStartAligned(nAnimals>=3,:,:),3);
dpMean      = nanmean(allDpStartAligned(nAnimals>=3,:,:),3);
dpSem       = sem(allDpStartAligned(nAnimals>=3,:,:),3); 
rxnAvgMean  = nanmean(allRxnAvgStartAligned(nAnimals>=3,:,:),3);
rxnAvgSem   = sem(allRxnAvgStartAligned(nAnimals>=3,:,:),3); 
rxnSdMean   = nanmean(allRxnSdStartAligned(nAnimals>=3,:,:),3);
rxnSdSem    = sem(allRxnSdStartAligned(nAnimals>=3,:,:),3); 
x           = repmat([1:maxSessions]',1,num_trial_types);

first200 = squeeze(allDpStartAligned(1,2,:));
[last200,first100,last100] = deal(nan*first200);
[firstLastSame200,firstLastSame100] = deal(false(size(allDpStartAligned,3),1));
for i = 1:size(allDpStartAligned,3)
    last200Idx = find(~isnan(allDpStartAligned(:,2,i)),1,'Last');
    if ~isempty(last200Idx)
        last200(i) = allDpStartAligned(last200Idx,2,i);
    end
    if last200Idx == 1
        firstLastSame200(i) = true;
    end
    
    first100Idx = find(~isnan(allDpStartAligned(:,3,i)),1,'First');
    if ~isempty(first100Idx)
        first100(i) = allDpStartAligned(first100Idx,3,i);
    end
    
    last100Idx = find(~isnan(allDpStartAligned(:,3,i)),1,'Last');
    if last100Idx == first100Idx
        firstLastSame100(i) = true;
    end
    if ~isempty(last100Idx)
        last100(i) = allDpStartAligned(last100Idx,3,i);
    end
    
end

% Save out if desired
% save([backdir(which('DalgleishHausser2020_All2PLearning_raw.mat'),1) filesep 'DalgleishHausser2020_All2PLearning_analysed.mat'],...
%     'lastPr','lastRxnMean','lastRxnSd','prMean','prSem','dpMean','dpSem',...
%     'first200','last200','first100','last100','firstLastSame100',...
%     'firstLastSame200','rxnAvgMean','rxnAvgSem','rxnSdMean','rxnSdSem',...
%     'maxSessions','x','animal_ids','all_pr','all_rxn_avg','all_rxn_sd',...
%     'all_dp','bhv_eg')

%% Figure 1 - figure supplement 3
% 1P training paradigm

load('DalgleishHausser2020_All1PLearning_raw.mat')

% Find stim pr, dp, rxn for each power
[aw_pr,aw_dp,aw_rxn_avg,aw_rxn_sd,aw_days_powers] = deal(cell(numel(all_bhv),1));
num_animals = numel(all_bhv);

% staircase sessions, organised by animals id/session id
for a = 1:num_animals
    for d = 1:numel(all_bhv(a).bhv_daily)
        ng = pyb_return_trials(all_bhv(a).bhv_daily{d},6,0);
        all_stim_vars = all_bhv(a).stimvars_daily{d};
        
        for s = 1:size(all_stim_vars,2)
            o = pyb_return_trials(all_bhv(a).bhv_daily{d},all_stim_vars(1,s),all_stim_vars(2,s));
            
            powers = [o(:).power];
            uniquePowers = sort(unique(powers),'Descend');
            
            for p = 1:numel(uniquePowers)
                theseTrials = powers==uniquePowers(p);
                thisPower = uniquePowers(p);
                if thisPower>10
                    aw_days_powers{a} = [aw_days_powers{a} ; d nan];
                else
                    aw_days_powers{a} = [aw_days_powers{a} ; d thisPower];
                end
                aw_pr{a} = [aw_pr{a} ; mean([o(theseTrials).responded])];
                aw_dp{a} = [aw_dp{a} ; dprime(sum([o(theseTrials).responded]),sum(theseTrials),sum([ng(:).responded]),numel(ng))];
                if sum([o(theseTrials).responded]) >= 3
                    aw_rxn_avg{a} = [aw_rxn_avg{a} ; nanmean([o(theseTrials).rxntime])];
                    aw_rxn_sd{a} = [aw_rxn_sd{a} ; nanstd([o(theseTrials).rxntime])];
                else
                    aw_rxn_avg{a} = [aw_rxn_avg{a} ; nan];
                    aw_rxn_sd{a} = [aw_rxn_sd{a} ; nan];
                end
            end
        end
    end
    flag = isnan(aw_days_powers{a}(:,2));
    aw_days_powers{a}(flag,:) = [];
    aw_pr{a}(flag) = [];
    aw_dp{a}(flag) = [];
    aw_rxn_avg{a}(flag) = [];
    aw_rxn_sd{a}(flag) = [];
end

% find session-wise power and dp
all_stims = [5 6];
[sw_pr,sw_mean_dp,sw_mean_power,sw_min_dp,sw_min_power] = deal(cell(num_animals,1));
for a = 1:num_animals
    hpSess = 0;
    for d = 1:numel(all_bhv(a).bhv_daily)
        ng = pyb_return_trials(all_bhv(a).bhv_daily{d},6,0);
        o = pyb_return_trials(all_bhv(a).bhv_daily{d},5,[]);
        if all(~isnan(unique([o(:).power])))
            hpSess = hpSess+1;
            for s = 1:size(all_stims,2)
                o = pyb_return_trials(all_bhv(a).bhv_daily{d},all_stims(s),[]);
                thesePowers = [o(:).power];
                sw_pr{a}(hpSess,s) = mean([o(:).responded]);
                sw_mean_power{a}(hpSess,s) = nanmean(thesePowers);
                sw_mean_dp{a}(hpSess,s) = dprime(sum([o(:).responded]),numel(o),sum([ng(:).responded]),numel(ng));
                
                uniquePowers = unique(thesePowers);
                nTrialsPowers = histc(thesePowers,uniquePowers);
                minPower = nanmin(uniquePowers(nTrialsPowers>20));
                sw_min_power{a}(hpSess,s) = minPower;
                theseTrials = thesePowers == sw_min_power{a}(hpSess,s);
                sw_min_dp{a}(hpSess,s) = dprime(sum([o(theseTrials).responded]),numel(o),sum([ng(:).responded]),numel(ng));
            end
        end
    end
end

maxNsessions = max(cellfun(@(x) size(x,1),sw_pr));
[sw_mean_dp_cat,sw_mean_power_cat,sw_min_dp_cat,sw_min_power_cat] = deal(nan(num_animals,maxNsessions));
for a = 1:num_animals
    sw_mean_dp_cat(a,1:size(sw_mean_dp{a},1)) = sw_mean_dp{a}(:,1);
    sw_mean_power_cat(a,1:size(sw_mean_power{a},1)) = sw_mean_power{a}(:,1);
    
    sw_min_dp_cat(a,1:size(sw_min_dp{a},1)) = sw_min_dp{a}(:,1);
    sw_min_power_cat(a,1:size(sw_min_power{a},1)) = sw_min_power{a}(:,1);
end

% calculate sliding performance and rxn
windowSize      = 20;
windowSizeCatch = 20;
[catStimRxn,catStimRxnCatch,catStimR,catStimRCatch] = deal(cell(num_animals,1));
for a = 1:num_animals
    for d = 1:numel(all_bhv(a).bhv_daily)
        o = pyb_return_trials(all_bhv(a).bhv_daily{d},5,[]);
        catStimR{a} = [catStimR{a} ; [o(:).responded]'];
        catStimRxn{a} = [catStimRxn{a} ; [o(:).rxntime]'];
        
        o = pyb_return_trials(all_bhv(a).bhv_daily{d},6,0);
        catStimRCatch{a} = [catStimRCatch{a} ; [o(:).responded]'];
        catStimRxnCatch{a} = [catStimRxnCatch{a} ; [o(:).rxntime]'];
    end
end
lens = cellfun(@(x) numel(x),catStimR);
maxLen = max(lens);
lensCatch = cellfun(@(x) numel(x),catStimRCatch);
maxLenCatch = max(lensCatch);
[slidingPr,slidingRxn,slidingRxnSd] = deal(nan(num_animals,maxLen));
[slidingPrCatch,slidingRxnCatch,slidingRxnSdCatch] = deal(nan(num_animals,maxLenCatch));
for a = 1:num_animals
    slidingPr(a,1:lens(a)) = movmean(catStimR{a},windowSize);
    slidingPrCatch(a,1:lensCatch(a)) = movmean(catStimRCatch{a},windowSize);
    
    slidingRxn(a,1:lens(a)) = movmean(catStimRxn{a},windowSize,'omitnan');
    slidingRxnCatch(a,1:lensCatch(a)) = movmean(catStimRxnCatch{a},windowSizeCatch,'omitnan');
    
    slidingRxnSd(a,1:lens(a)) = movstd(catStimRxn{a},windowSize,'omitnan');
    slidingRxnSdCatch(a,1:lensCatch(a)) = movstd(catStimRxnCatch{a},windowSizeCatch,'omitnan');
end

% get performance per day per power
allUniquePowers = unique(cell2mat(cellfun(@(x) unique(x(:,2)),aw_days_powers,'UniformOutput',0)));
allUniquePowers = sort(allUniquePowers(allUniquePowers>0),'Descend'); %[10 7.5 5 2.5 1.75 1.25 1 0.5];
firstPowers = cell2mat(cellfun(@(x) x(1,:),aw_days_powers,'UniformOutput',0));
[all_days,all_pr,all_dp] = deal(zeros(num_animals,numel(allUniquePowers)));
for p = 1:numel(allUniquePowers)
    for a = 1:num_animals
        theseIdcs = aw_days_powers{a}(:,2) == allUniquePowers(p);
        all_days(a,p) = mean(aw_days_powers{a}(theseIdcs,1));
        all_pr(a,p) = mean(aw_pr{a}(theseIdcs,1));
        all_dp(a,p) = mean(aw_dp{a}(theseIdcs,1));
    end
end

% Analyse high/low power psychometrics
highPowers  = [0.25 0.2 0.15 0.1 0.05 0];
lowPowers   = [0.1 0.08 0.06 0.04 0.02 0];
nHighPowers = numel(highPowers);
nLowPowers  = numel(lowPowers);
[maxSessionsHighPower,maxSessionsLowPower] = deal(0);
for a = 1:numel(psych_bhv)
    psych_bhv(a).highPowerSessions = [];
    psych_bhv(a).lowPowerSessions = [];
    for d = 1:numel(psych_bhv(a).powers)
        thesePowers = psych_bhv(a).powers{d}(~isnan(psych_bhv(a).powers{d}));
        if all(ismember(highPowers,thesePowers) & ismember(thesePowers,highPowers))
            psych_bhv(a).highPowerSessions = [psych_bhv(a).highPowerSessions ; d];
        elseif all(ismember(lowPowers,thesePowers) & ismember(thesePowers,lowPowers))
            psych_bhv(a).lowPowerSessions = [psych_bhv(a).lowPowerSessions ; d];
        end
    end
    maxSessionsHighPower = max([maxSessionsHighPower numel(psych_bhv(a).highPowerSessions)]);
    maxSessionsLowPower = max([maxSessionsLowPower numel(psych_bhv(a).lowPowerSessions)]);
end

[crossDayHighPsychsPr,crossDayHighPsychN,crossDayHighPsychR,crossDayHighPsychPower,...
    crossDayHighPsychsRxnMean,crossDayHighPsychsRxnSd] = ...
    deal(nan(numel(psych_bhv),nHighPowers,maxSessionsHighPower));
[crossDayLowPsychsPr,crossDayLowPsychN,crossDayLowPsychR,crossDayLowPsychPower,...
    crossDayLowPsychsRxnMean,crossDayLowPsychsRxnSd] = ...
    deal(nan(numel(psych_bhv),nHighPowers,maxSessionsLowPower));
for a = 1:numel(psych_bhv)
    % high power
    hpSess = psych_bhv(a).highPowerSessions;
    for d = 1:numel(hpSess)
        for p = 1:numel(highPowers)
            theseTrials = [psych_bhv(a).bhv_daily{hpSess(d)}(:).power] == highPowers(p);
            crossDayHighPsychN(a,p,d) = sum(theseTrials);
            crossDayHighPsychR(a,p,d) = sum([psych_bhv(a).bhv_daily{hpSess(d)}(theseTrials).responded]);
            crossDayHighPsychPower(a,p,d) = highPowers(p);
            crossDayHighPsychsPr(a,p,d) = mean([psych_bhv(a).bhv_daily{hpSess(d)}(theseTrials).responded]);
            
            if crossDayHighPsychR(a,p,d)>=3
                crossDayHighPsychsRxnMean(a,p,d) = nanmean([psych_bhv(a).bhv_daily{hpSess(d)}(theseTrials).rxntime]);
                crossDayHighPsychsRxnSd(a,p,d) = nanstd([psych_bhv(a).bhv_daily{hpSess(d)}(theseTrials).rxntime]);
            end
        end
    end
    
    lpSess = psych_bhv(a).lowPowerSessions;
    for d = 1:numel(lpSess)
        for p = 1:numel(lowPowers)
            theseTrials = [psych_bhv(a).bhv_daily{lpSess(d)}(:).power] == lowPowers(p);
            crossDayLowPsychN(a,p,d) = sum(theseTrials);
            crossDayLowPsychR(a,p,d) = sum([psych_bhv(a).bhv_daily{lpSess(d)}(theseTrials).responded]);
            crossDayLowPsychPower(a,p,d) = lowPowers(p);
            crossDayLowPsychsPr(a,p,d) = mean([psych_bhv(a).bhv_daily{lpSess(d)}(theseTrials).responded]);
            
            if crossDayLowPsychN(a,p,d)>=3
                crossDayLowPsychsRxnMean(a,p,d) = nanmean([psych_bhv(a).bhv_daily{lpSess(d)}(theseTrials).rxntime]);
                crossDayLowPsychsRxnSd(a,p,d) = nanstd([psych_bhv(a).bhv_daily{lpSess(d)}(theseTrials).rxntime]);
            end
        end
    end
end

% high power first and last session
tmp = squeeze(crossDayHighPsychsPr(:,highPowers==0.05,:));
firstLastHighPsych = nan(numel(psych_bhv),2);
for a = 1:numel(psych_bhv)
    idcs = find(~isnan(tmp(a,:)));
    if numel(idcs)>1
        firstLastHighPsych(a,:) = tmp(a,[idcs(1) idcs(end)]);
    end
end
firstLastHighPsych(any(isnan(firstLastHighPsych),2),:) = [];

% do stats on low power (avg of first 3 sessions)
lowPowerPsychAvg = squeeze(nanmean(crossDayLowPsychsPr(:,:,1:3),3));
lowPowerPsychAvgRxnMean = squeeze(nanmean(crossDayLowPsychsRxnMean(:,:,1:3),3));
lowPowerPsychAvgRxnSd = squeeze(nanmean(crossDayLowPsychsRxnSd(:,:,1:3),3));

lowPowerPsychAvg(any(isnan(lowPowerPsychAvg),2),:) = [];
lowPowerPsychAvgRxnMean(all(isnan(lowPowerPsychAvgRxnMean),2),:) = [];
lowPowerPsychAvgRxnSd(all(isnan(lowPowerPsychAvgRxnSd),2),:) = [];
pr_pVal = size(lowPowerPsychAvg,2)-1;
for i = 1:size(lowPowerPsychAvg ,2)-1
    [pr_pVal(i),~] = signrank(lowPowerPsychAvg(:,end),lowPowerPsychAvg(:,i));
end
pr_pVal = pr_pVal * size(pr_pVal,2);
pr_sig = pr_pVal < 0.05;

lowPowerPrSt = qstat(lowPowerPsychAvg(:,[find(lowPowers==0.1) find(lowPowers==0.02)]));
lowPowerRxnMeanSt = qstat(lowPowerPsychAvgRxnMean(:,[find(lowPowers==0.1) find(lowPowers==0.02)]));
lowPowerRxnSdSt = qstat(lowPowerPsychAvgRxnSd(:,[find(lowPowers==0.1) find(lowPowers==0.02)]));

% find no. 1P sessions of each type
nSessionsStaircase = sum(~isnan(sw_mean_dp_cat),2);
nSessionsHighPsych = sum(squeeze(all(~isnan(crossDayHighPsychsPr),2)),2);
nSessionsLowPsych = sum(squeeze(all(~isnan(crossDayLowPsychsPr),2)),2);

% Save out if desired
% save([backdir(which('DalgleishHausser2020_All2PLearning_raw.mat'),1) filesep 'DalgleishHausser2020_All1PLearning_analysed.mat',...
%     'nSessionsStaircase','nSessionsHighPsych','nSessionsLowPsych',...
%     'slidingPr','slidingPrCatch','slidingRxnCatch','slidingRxn',...
%     'slidingRxnSdCatch','slidingRxnSd','sw_mean_dp_cat',...
%     'sw_mean_power_cat','sw_min_dp_cat','sw_min_power_cat',...
%     'maxSessionsHighPower','crossDayHighPsychsPr','highPowers',...
%     'crossDayLowPsychsPr','lowPowers','crossDayLowPsychsRxnMean',...
%     'crossDayLowPsychsRxnSd','animal_ids')

