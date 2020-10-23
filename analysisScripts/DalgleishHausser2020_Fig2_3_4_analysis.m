%% Process data (Figures 2 - 4)

%%%%%%%%% Process (necessary for all figures 2 - 4)
fprintf(['- Loading data...\n'])
if ~exist('tw_response_zs') || ~exist('tw_traces_zs')
    fprintf(['- Importing saved data...\n'])
    load('DalgleishHausser2020_imaging_raw.mat')
else
    fprintf(['- Using data in workspace...\n'])
end

fprintf(['- Analysing target and network response metrics...\n'])
rng('default')
metric                  = 'proportion';     % proportion | amplitude
num_shuffles            = 100;              % shuffles for hit/miss matching
max_var                 = num_cells==200;
ng_var                  = num_cells==0;
thold_var               = num_cells==50;

% quantify target/network response
[e_all_matched_n,i_all_matched_n,e_off_matched_n,i_off_matched_n,e_t_matched_n,i_t_matched_n,...
    e_all_n,e_all_hit_n,e_all_miss_n,i_all_n,i_all_hit_n,i_all_miss_n,...
    e_off_n,e_off_hit_n,e_off_miss_n,i_off_n,i_off_hit_n,i_off_miss_n,...
    e_t_n,e_t_hit_n,e_t_miss_n,i_t_n,i_t_hit_n,i_t_miss_n,...
    matched_num_trials,p_targets,...
    all_e_t_n,all_e_t_n_catch,...
    e_t_n_catch,e_t_n_hit_catch,e_t_n_miss_catch,...
    i_t_n_catch,i_t_n_hit_catch,i_t_n_miss_catch,...
    e_t_matched_n_catch,i_t_matched_n_catch] = ...
    deal(zeros(num_sessions,num_stim_vars));
pop_size = zeros(num_sessions,1);
[cw_e_all_matched,cw_i_all_matched,targ_z_dist,targeted_expressing_flag,...
    targets_with_neurons,untargeted_expressing_flag,untargeted_nonexpressing_flag] ...
    = deal(cell(num_sessions,num_stim_vars));
for s = 1:num_sessions
    
    % these targets and pop size
    pop_size(s,1)                   = size(tw_response_zs{s,1},1);

    for t = 1:num_stim_vars

        % return hit/miss idcs for stim and catch trials
        hit_idcs                    = find(tw_response{s,t});
        n_hits                      = numel(hit_idcs);
        miss_idcs                   = find(~tw_response{s,t});
        n_misses                    = numel(miss_idcs);
        
        % STIM TRIALS   
        % all
        tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric);
        e_all_n(s,t)                = mean(nansum(tmp,1));
        e_all_hit_n(s,t)            = mean(nansum(tmp(:,hit_idcs),1));
        e_all_miss_n(s,t)           = mean(nansum(tmp(:,miss_idcs),1));

        tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric);
        i_all_n(s,t)                = mean(nansum(tmp,1));
        i_all_hit_n(s,t)            = mean(nansum(tmp(:,hit_idcs),1));
        i_all_miss_n(s,t)           = mean(nansum(tmp(:,miss_idcs),1));
        
        % targets
        this_flag                   = ~targeted_flag{s,t};
        tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric,'filter',this_flag);
        e_t_n(s,t)                  = mean(nansum(tmp,1));
        e_t_hit_n(s,t)              = mean(nansum(tmp(:,hit_idcs),1));
        e_t_miss_n(s,t)             = mean(nansum(tmp(:,miss_idcs),1));
        
        tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric,'filter',this_flag);
        i_t_n(s,t)                  = mean(nansum(tmp,1));
        i_t_hit_n(s,t)              = mean(nansum(tmp(:,hit_idcs),1));
        i_t_miss_n(s,t)             = mean(nansum(tmp(:,miss_idcs),1));
        
        all_e_t_n(s,t)              = mean(sum(tw_response_zs{s,t}(targeted_flag{s,t},:) > pos_thold{s}(targeted_flag{s,t}),1),2);
        all_e_t_n_catch(s,t)        = mean(sum(tw_response_zs{s,ng_var}(targeted_flag{s,t},:) > pos_thold{s}(targeted_flag{s,t}),1),2);  
             
        % target catch rates
        if t == find(ng_var)
            for sv = 1:num_stim_vars
                this_flag                   = ~targeted_flag{s,sv};
                tmp                         = proc_tw_response(tw_response_zs{s,ng_var},pos_thold{s},'above',metric,'filter',this_flag);
                e_t_n_catch(s,sv)           = mean(nansum(tmp,1));
                e_t_n_hit_catch(s,sv)       = mean(nansum(tmp(:,hit_idcs),1));
                e_t_n_miss_catch(s,sv)      = mean(nansum(tmp(:,miss_idcs),1));
                
                tmp                         = proc_tw_response(tw_response_zs{s,ng_var},neg_thold{s},'below',metric,'filter',this_flag);
                i_t_n_catch(s,sv)           = mean(nansum(tmp,1));
                i_t_n_hit_catch(s,sv)       = mean(nansum(tmp(:,hit_idcs),1));
                i_t_n_miss_catch(s,sv)      = mean(nansum(tmp(:,miss_idcs),1));
            end
        end
        
        % background
        this_flag                   = targeted_flag{s,max_var};
        tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric,'filter',this_flag);
        e_off_n(s,t)                = mean(nansum(tmp,1));
        e_off_hit_n(s,t)            = mean(nansum(tmp(:,hit_idcs),1));
        e_off_miss_n(s,t)           = mean(nansum(tmp(:,miss_idcs),1));
        
        tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric,'filter',this_flag);
        i_off_n(s,t)                = mean(nansum(tmp,1));
        i_off_hit_n(s,t)            = mean(nansum(tmp(:,hit_idcs),1));
        i_off_miss_n(s,t)           = mean(nansum(tmp(:,miss_idcs),1));
        
        % shuffle to match hit/miss numbers
        [cw_e_all_matched{s,t},cw_i_all_matched{s,t}] = deal(0*tw_response_zs{s,t}(:,1));
        for i = 1:num_shuffles
            idcs2use = [];
            % stim trials
            % NB in condition with only 1 response type (only hits, only
            % misses) this will return [] because randperm contingent on no
            % trials of deficient type (i.e. 0 hits if only misses)
            if (n_hits > n_misses) % permute hits
                idcs2use = [miss_idcs ; hit_idcs(randperm(n_hits,n_misses))];
            elseif (n_hits < n_misses) % permute misses
                idcs2use = [hit_idcs ; miss_idcs(randperm(n_misses,n_hits))];
            elseif n_hits == n_misses % just use all trials (NB this used to just be ELSE)
                idcs2use = [hit_idcs ; miss_idcs];
            end
            
            % STIM TRIALS    
            % all
            tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric);
            e_all_matched_n(s,t)        = e_all_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric);
            i_all_matched_n(s,t)        = i_all_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            % targets
            this_flag                   = ~targeted_flag{s,t};
            tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric,'filter',this_flag);
            e_t_matched_n(s,t)          = e_t_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric,'filter',this_flag);
            i_t_matched_n(s,t)          = i_t_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            if t == find(ng_var)
                for sv = 1:num_stim_vars
                    this_flag                   = ~targeted_flag{s,sv};
                    tmp                         = proc_tw_response(tw_response_zs{s,ng_var},pos_thold{s},'above',metric,'filter',this_flag);
                    e_t_matched_n_catch(s,sv)   = e_t_matched_n_catch(s,sv) + mean(nansum(tmp(:,idcs2use),1));
                    
                    tmp                         = proc_tw_response(tw_response_zs{s,ng_var},neg_thold{s},'below',metric,'filter',this_flag);
                    i_t_matched_n_catch(s,sv)   = i_t_matched_n_catch(s,sv) + mean(nansum(tmp(:,idcs2use),1));

                end
            end
            
            % background
            this_flag                   = targeted_flag{s,max_var};
            tmp                         = proc_tw_response(tw_response_zs{s,t},pos_thold{s},'above',metric,'filter',this_flag);
            e_off_matched_n(s,t)        = e_off_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            tmp                         = proc_tw_response(tw_response_zs{s,t},neg_thold{s},'below',metric,'filter',this_flag);
            i_off_matched_n(s,t)        = i_off_matched_n(s,t) + mean(nansum(tmp(:,idcs2use),1));
            
            % cell-wise matched responses
            cw_e_all_matched{s,t}       = cw_e_all_matched{s,t} + mean(tw_response_zs{s,t}(:,idcs2use)>pos_thold{s},2);
            cw_i_all_matched{s,t}       = cw_i_all_matched{s,t} + mean(tw_response_zs{s,t}(:,idcs2use)<neg_thold{s},2);
            
            if mod(i,500)==0
                fprintf(['Session ' num2str(s) ', perm ' num2str(i) '\n'])
            end
               
        end
        cw_e_all_matched{s,t}           = cw_e_all_matched{s,t} / num_shuffles;
        cw_i_all_matched{s,t}           = cw_i_all_matched{s,t} / num_shuffles;
        
        matched_num_trials(s,t)         = numel(idcs2use);
    end
end
mean_num_trials = mean(matched_num_trials,1);
pr = cellfun(@(x) mean(x),tw_response);
allTargets = cellfun(@(x) sum(x),targeted_flag);
targeted = repmat(num_cells,size(allTargets,1),1);

% average across shuffles
e_all_matched_n             = e_all_matched_n / num_shuffles;
i_all_matched_n             = i_all_matched_n / num_shuffles;
e_t_matched_n               = e_t_matched_n / num_shuffles;
i_t_matched_n               = i_t_matched_n / num_shuffles;
e_t_matched_n_catch         = e_t_matched_n_catch / num_shuffles;
i_t_matched_n_catch         = i_t_matched_n_catch / num_shuffles;
e_off_matched_n             = e_off_matched_n / num_shuffles;
i_off_matched_n             = i_off_matched_n / num_shuffles;

% convert to proportions
all_e_t                     = all_e_t_n./pop_size;
all_e_t_catch               = all_e_t_n_catch./pop_size;

e_all                       = e_all_n./pop_size;
e_all_hit                   = e_all_hit_n./pop_size;
e_all_miss                  = e_all_miss_n./pop_size;
i_all                       = i_all_n./pop_size;
i_all_hit                   = i_all_hit_n./pop_size;
i_all_miss                  = i_all_miss_n./pop_size;
e_t                         = e_t_n./pop_size;
e_t_hit                     = e_t_hit_n./pop_size;
e_t_miss                    = e_t_miss_n./pop_size;
e_t_catch                   = e_t_n_catch./pop_size;
e_t_hit_catch               = e_t_n_hit_catch./pop_size;
e_t_miss_catch              = e_t_n_miss_catch./pop_size;
i_t                         = i_t_n./pop_size;
i_t_hit                     = i_t_hit_n./pop_size;
i_t_miss                    = i_t_miss_n./pop_size;
e_off                       = e_off_n./pop_size;
e_off_hit                   = e_off_hit_n./pop_size;
e_off_miss                  = e_off_miss_n./pop_size;
i_off                       = i_off_n./pop_size;
i_off_hit                   = i_off_hit_n./pop_size;
i_off_miss                  = i_off_miss_n./pop_size;

e_all_matched               = e_all_matched_n./pop_size;
i_all_matched               = i_all_matched_n./pop_size;
e_t_matched                 = e_t_matched_n./pop_size;
i_t_matched                 = i_t_matched_n./pop_size;
e_t_matched_catch           = e_t_matched_n_catch./pop_size;
i_t_matched_catch           = i_t_matched_n_catch./pop_size;
e_off_matched               = e_off_matched_n./pop_size;
i_off_matched               = i_off_matched_n./pop_size;

%%%%%%%%% Pyschometric curve analyses (Figure 2)
fprintf(['- Fitting psychometric curves...\n'])
% get data
nCellsPsych = e_t_n(:,~ng_var);
nResponses  = cellfun(@(x) sum(x),tw_response(:,~ng_var));
nTrials     = cellfun(@(x) numel(x),tw_response(:,~ng_var));
pResponse   = nResponses./nTrials;
xScatter    = nCellsPsych;
pts2analyse = [0.1 0.5 0.9];

% remove problematic data
flag = isnan(nCellsPsych) | isinf(nCellsPsych) | nCellsPsych==0 | nCellsPsych<1;
nCellsPsych(flag) = nan;

% prepare structure (x | nCorrect | total) options
data                = [nCellsPsych(:) nResponses(:) nTrials(:)];
data(any(isnan(data),2),:) = [];
guessRate           = mean(pr(:,ng_var));   % global
lapseRate           = mean(max(pr(:,~ng_var),[],2));  % global
options             = struct;
options.sigmoidName = 'logn'; 
options.expType     = 'YesNo';
options.expN        = 1;

% GLOBAL FIT
options.fixedPars   = [nan ; nan ; 1-lapseRate ; guessRate ; nan];
numCellsBhvFit      = psignifit(data,options);
numCellsBhvFitSlope = getSlopePC(numCellsBhvFit,0.5);
numCellsBhvFitR2    = calculate_r2(data(:,2)./data(:,3),psigniGetCurve(numCellsBhvFit,data(:,1)));
numCellsBhvFitPoints = 0*pts2analyse;
numCellsBhvFitPointsCI = zeros(numel(pts2analyse),2);
for i = 1:numel(pts2analyse)
    options.threshPC = pts2analyse(i);
    tmp = psignifit(data,options);
    numCellsBhvFitPoints(i) = exp(tmp.Fit(1));
    numCellsBhvFitPointsCI(i,:) = exp(tmp.conf_Intervals(1,:,1));
end

% INDIVIDUAL FITS
num_sessions = size(e_t_n,1);
options.threshPC = 0.5;
clear numCellsBhvFitInd;
xFit = linspace(min(nCellsPsych(:)),max(nCellsPsych(:)),1000);
psycCurves = zeros(num_sessions,length(xFit));
[numCellsBhvFitIndR2,numCellsBhvIndSlopes] = deal(zeros(num_sessions,1));
numCellsBhvFitIndPoints = zeros(num_sessions,numel(pts2analyse));
for s = 1:num_sessions
    
    data            = [nCellsPsych(s,:)' nResponses(s,:)' nTrials(s,:)'];
    data(any(isnan(data),2),:) = [];
    options.fixedPars = ...
                    [nan ; nan ; 1-max(pr(s,~ng_var)) ; pr(s,ng_var) ; nan];
    numCellsBhvFitInd(s) = psignifit(data,options);
    psycCurves(s,:) = psigniGetCurve(numCellsBhvFitInd(s),xFit);
    numCellsBhvFitIndR2(s) = calculate_r2(data(:,2)./data(:,3),psigniGetCurve(numCellsBhvFitInd(s),data(:,1)));
    for i = 1:numel(pts2analyse)
        numCellsBhvFitIndPoints(s,i) = getThreshold(numCellsBhvFitInd(s),pts2analyse(i),true);
    end
    
    numCellsBhvIndSlopes(s) = getSlopePC(numCellsBhvFitInd(s),0.5);
end

%%%%%%%%% Distance vs activation (Figure 3)
fprintf(['- Analysing distance vs activation...\n'])
[cw_e_off_matched,cw_i_off_matched,nearestTarget_off,...
    nearestTarget_ed,nearestTarget_xy,nearestTarget_z,...
    cw_e_all_matched_sub,cw_i_all_matched_sub,cw_e_t_matched_sub,...
    cw_e_off_matched_sub,cw_i_off_matched_sub] ...
    = deal(cell(num_sessions,num_stim_vars));
for s = 1:num_sessions
    for t = 1:num_stim_vars
        % split targets and network
        cw_e_off_matched{s,t}           = cw_e_all_matched{s,t}(~targeted_flag{s,t});
        cw_i_off_matched{s,t}           = cw_i_all_matched{s,t}(~targeted_flag{s,t});
        
        % find nearest stimulated cell/target
        if num_cells(t)>0
            % find euclean, xy and z distances of neurons to targets
            d                       = euclidean_distance(session(s).s2p.centroids_proc,session(s).targets(t).targets_xyz_proc);
            xyd                     = euclidean_distance(session(s).s2p.centroids_proc(:,1:2),session(s).targets(t).targets_xyz_proc(:,1:2));
            zd                      = round(session(s).s2p.centroids_proc(:,end) - session(s).targets(t).targets_xyz_proc(:,end)',1);
            
            % find closest euclidean target
            nearestTarget_off{s,t} = min(d(~targeted_flag{s,t},:),[],2);
            
            % find xy and z distances from this closest target
            [nearestTarget_ed{s,t},minI]   ...
                = min(d,[],2);
            nearestTarget_xy{s,t}   = diag(xyd(:,minI));
            nearestTarget_z{s,t}    = diag(zd(:,minI));
        else
            nearestTarget_off{s,t}  = [];
        end
        
        % subtract off spontaneous rates
        if t ~= find(ng_var)
            cw_e_all_matched_sub{s,t}   = cw_e_all_matched{s,t} - cw_e_all_matched{s,ng_var};
            cw_i_all_matched_sub{s,t}   = cw_i_all_matched{s,t} - cw_i_all_matched{s,ng_var};
            
            cw_e_off_matched_sub{s,t}   = cw_e_all_matched_sub{s,t}(~targeted_flag{s,t});
            cw_i_off_matched_sub{s,t}   = cw_i_all_matched_sub{s,t}(~targeted_flag{s,t});
        end
    end
end

% xy vs p(response) - avg across animals
xyBins      = [0:10:110];
absDepth    = true;
zDepths     = [0 33.3 66.6 99.9];

% lateral distance vs activation (collapsed across axial)
% [trial_type lateral_bin animal]
[xyE,xyI,xyECatch,xyICatch] = deal(nan(num_stim_vars-1,numel(xyBins)-1,num_sessions));
for s = 1:num_sessions
    for t = 2:num_stim_vars
        binIdcs = discretize(nearestTarget_xy{s,t},xyBins);
        binIdcs(isnan(binIdcs)) = 0;
        
        for j = 1:numel(xyBins)-1
            theseIdcs = binIdcs == j;
            xyE(t-1,j,s) = nanmean(cw_e_all_matched_sub{s,t}(theseIdcs));
            xyI(t-1,j,s) = nanmean(cw_i_all_matched_sub{s,t}(theseIdcs));
            
            xyECatch(t-1,j,s) = nanmean(cw_e_all_matched{s,ng_var}(theseIdcs));
            xyICatch(t-1,j,s) = nanmean(cw_i_all_matched{s,ng_var}(theseIdcs));
        end
    end
end

% axial & lateral distance vs activation
if ~absDepth
    zDepths = unique(sort([zDepths -zDepths],'Descend'));
end
[xyzE,xyzI] = deal(nan(numel(zDepths),numel(xyBins)-1,num_sessions));
[tmpE,tmpI] = deal(nan(1,num_stim_vars-1));
for s = 1:num_sessions
    for i = 1:numel(zDepths)
        for j = 1:numel(xyBins)-1
            for t = 2:num_stim_vars
                binIdcs = discretize(nearestTarget_xy{s,t},xyBins);
                binIdcs(isnan(binIdcs)) = 0;
                if absDepth
                    theseIdcs = (binIdcs == j) & (abs(nearestTarget_z{s,t}) == zDepths(i));
                else
                    theseIdcs = (binIdcs == j) & (nearestTarget_z{s,t} == zDepths(i));
                end
                tmpE(t-1) = nanmean(cw_e_all_matched_sub{s,t}(theseIdcs));
                tmpI(t-1) = nanmean(cw_i_all_matched_sub{s,t}(theseIdcs));
            end
            xyzE(i,j,s) = nanmean(tmpE);
            xyzI(i,j,s) = nanmean(tmpI);
        end
    end
end
xyBins = xyBins(1:end-1);

%%%%%%%%% Network response models, permuted (Figure 4)
fprintf(['- Fitting network response vs behaviour models...\n'])
% get data
xAll                = [];
fewFlag             = reshape(e_t_matched_n(:,~ng_var),[],1)<1;
xAll(:,end+1)       = reshape(e_t_matched_n(:,~ng_var),[],1);
xAll(:,end+1)       = reshape(e_off_matched(:,~ng_var),[],1);
xAll(:,end+1)       = reshape(i_off_matched(:,~ng_var),[],1);
nResponsesMult      = reshape(cellfun(@(x) sum(x),tw_response(:,~ng_var)),[],1);
nTrialsMult         = reshape(cellfun(@(x) numel(x),tw_response(:,~ng_var)),[],1);

% remove problematic data
flag = any(isnan(xAll) | isinf(xAll),2) | fewFlag;
xAll(flag,:) = []; nResponsesMult(flag) = []; nTrialsMult(flag) = [];

% fit options
lapseRate           = 1-mean(max(pr(:,~ng_var),[],2));
guessRate           = mean(pr(:,ng_var));
options             = struct;
options.sigmoidName = 'logn'; 
options.expType     = 'YesNo';
options.expN        = 1;
options.fixedPars   = [nan ; nan ; lapseRate ; guessRate ; nan];

% fit models (to all data, unpermuted)
clear 'predictorFits';
nPredictors = size(xAll,2);
r2 = zeros(1,nPredictors);
for i = 1:nPredictors
    data = [xAll(:,i) nResponsesMult nTrialsMult];
    predictorFits(i) = psignifit(data,options);
    yPred = psigniGetCurve(predictorFits(i),data(:,1));
    yAct = data(:,2)./data(:,3);  
    r2(i) = calculate_r2(yAct(:),yPred(:));
end

% Permute
try
    load('DalgleishHausser2020_Fig4_PermutedIndividualModel.mat')
catch
    clc
    spmd
        warning('off')
    end
    indModel_nP = 10000;
    n = size(xAll,1);
    rng('default')
    [indModel_r2PermFit,indModel_r2PermVal] = deal(zeros(indModel_nP,nPredictors));
    parfor i = 1:nPredictors
        data = [xAll(:,i) nResponsesMult nTrialsMult];
        for p = 1:indModel_nP
            [trainInd,valInd] = dividerand(n,0.8,0.2);
            tmp = psignifit(data(trainInd,:),options);
            indModel_r2PermFit(p,i) = calculate_r2(data(trainInd,2)./data(trainInd,3),psigniGetCurve(tmp,data(trainInd,1)));
            indModel_r2PermVal(p,i) = calculate_r2(data(valInd,2)./data(valInd,3),psigniGetCurve(tmp,data(valInd,1)));
            if mod(p,25)==0
                sprintf(['Predictor = ' num2str(i) ', Iter: ' num2str(p)])
            end
        end
    end
    spmd
        warning('on')
    end
    %save([backdir(which('DalgleishHausser2020_Fig4_PermutedIndividualModel.mat'),1) filesep 'DalgleishHausser2020_Fig4_PermutedIndividualModel.mat'],'indModel_r2PermFit','indModel_r2PermVal','indModel_nP')
end

isRun.figure2 = true;
isRun.figure3 = true;
isRun.figure4 = true;
