%% Load data
% NB change "base_dir" below to the directory containing raw data
% downloaded from: https://doi.org/10.6084/m9.figshare.13128950
%
% This data directory should contain the sessionDate_animalID directories
% unzipped from the above figshare database (see the "Raw data download and
% import" section of DalgleishHausser2020 readme):
%
% --- ~/DalgleishHausser2020Repo/rawDataDirectory/
%     |--- 20191203_L541/
%     |    --- experimentDataFile1
%     |    --- experimentDataFile2
%     |--- 20191009_L544/
%     |    --- experimentDataFile1
%     |    --- experimentDataFile2
%     . etc.
%
% It is recommended to use the unzipDirectories.m function included in this
% repo to ensure that the correct directory structure is retained on
% unzipping (see documentation for that function).

base_dir = '/Users/henrydalgleish/Dropbox/Analysis/CorrBhvExpt_LR/Imaging/NumCells/';

min_lick_t          = 0.15;
response_window     = 1.15;
easy_stim           = 7;
easy_var            = 1;
sated_thold         = 0.7;
min_dur             = 0.5;
sess_stop_window    = 10;

% find experiment directories
filter = '[0-9]{8}_[A-Z]{1,2}[0-9]{3}$';
[exp_dirs] = return_fullfile(base_dir,'re',filter);

% load all data
session = struct;
for d = 1:numel(exp_dirs)
    fprintf(['\nLoading ' exp_dirs{d} '\n'])
    session(d).dir = exp_dirs{d};
    paqs = sort_nat(return_fullfile(exp_dirs{d},'_paqanalysis.mat'));
    for p = 1:numel(paqs)
        tmp = load(paqs{p});
        session(d).paq(p) = paq_add_responses(tmp.pa,'response_window',response_window);
    end
    [session(d).s2p.f,session(d).s2p.np,session(d).s2p.centroids,session(d).s2p.stat,session(d).s2p.ops,session(d).s2p.spks] = load_s2py_stack(fullfile(exp_dirs{d},'Fall.mat'));
    session(d).targets = load_tpbs_targets(fullfile(exp_dirs{d},'targets'));
end

% Concatenate paq files
num_planes = numel(unique(session(1).s2p.centroids(:,end)));
for i = 1:numel(session)
    
    % process paq frames accoring to s2p
    offsets = cumsum([0 double(session(i).s2p.ops.frames_per_file(1:end-1))]);
    for p = 1:numel(session(i).paq)
        session(i).paq(p).cat_offset = offsets(p);
        for s = 1:numel(session(i).paq(p).stims)
            if ~isempty(session(i).paq(p).stims(s).in) && ~isempty(session(i).paq(p).stims(s).in.frames{1})
                session(i).paq(p).stims(s).in.frames_proc = cellfun(@(x) round(x/num_planes),session(i).paq(p).stims(s).in.frames,'UniformOutput',0);
                session(i).paq(p).stims(s).in.frames_offset = cellfun(@(x) x+session(i).paq(p).cat_offset,session(i).paq(p).stims(s).in.frames_proc,'UniformOutput',0);
            end
            if ~isempty(session(i).paq(p).stims(s).out) && ~isempty(session(i).paq(p).stims(s).out.frames{1})
                session(i).paq(p).stims(s).out.frames_proc = cellfun(@(x) round(x/num_planes),session(i).paq(p).stims(s).out.frames,'UniformOutput',0);
                session(i).paq(p).stims(s).out.frames_offset = cellfun(@(x) x+session(i).paq(p).cat_offset,session(i).paq(p).stims(s).out.frames_proc,'UniformOutput',0);
            end
        end
    end
    
    % concatenate stim frames, reward classification and reaction times
    stim_var_n_t_rxn_trial = [];
    for s = 1:numel(session(i).paq(p).stims)
        num_vars = numel(session(i).paq(1).stims(s).vars);
        if num_vars>0
            [session(i).cat_paq(s).stim_frames,session(i).cat_paq(s).responses,session(i).cat_paq(s).rxntimes,session(i).cat_paq(s).stim_frames_raw] = deal(cell(1,num_vars));
            for v = 1:num_vars
                for p = 1:numel(session(i).paq)
                    % use out triggers (i.e. 'FromPV')for stim trials
                    if ~isempty(session(i).paq(p).stims(s).out) && ~isempty(session(i).paq(p).stims(s).out.frames{v})
                        session(i).cat_paq(s).stim_frames{v} = [session(i).cat_paq(s).stim_frames{v} ; session(i).paq(p).stims(s).out.frames_offset{v}];
                    % else us in triggers (i.e. from PyB) for catch trials 
                    else
                        session(i).cat_paq(s).stim_frames{v} = [session(i).cat_paq(s).stim_frames{v} ; session(i).paq(p).stims(s).in.frames_offset{v}];
                    end
                    session(i).cat_paq(s).responses{v} = [session(i).cat_paq(s).responses{v} ; session(i).paq(p).stims(s).responses{v}'];
                    session(i).cat_paq(s).rxntimes{v} = [session(i).cat_paq(s).rxntimes{v} ; session(i).paq(p).stims(s).rxntimes{v}'];
                end
                stim_var_n_t_rxn_trial = [stim_var_n_t_rxn_trial ; 
                                          [[s v] .* ones(numel(session(i).cat_paq(s).stim_frames{v}),2) ...
                                          [1:numel(session(i).cat_paq(s).stim_frames{v})]' ...
                                           session(i).cat_paq(s).stim_frames{v} ...
                                           session(i).cat_paq(s).rxntimes{v}]];
            end
        end
    end
    % keep raw paqs to exclude all trials/trial licking (i.e. not filtered for early/demotivated trials)
    session(i).cat_paq_raw = session(i).cat_paq;
    
    [~,order] = sort(stim_var_n_t_rxn_trial(:,4),'Ascend');
    n_trials = size(stim_var_n_t_rxn_trial,1);
    session(i).stim_var_n_t_rxn_trial = [stim_var_n_t_rxn_trial(order,:) [1:n_trials]'];
    session(i).stim_var_n_t_rxn_trial_raw = session(i).stim_var_n_t_rxn_trial;
    
    % remove too quick trials and unmotivated trials at end of sessions
    easy_trials                 = session(i).stim_var_n_t_rxn_trial(session(i).stim_var_n_t_rxn_trial(:,1) == easy_stim & ...
                                                                    session(i).stim_var_n_t_rxn_trial(:,2) == easy_var,end)';
    first_unmotivated_trial     = find(movmean(session(i).cat_paq(easy_stim).responses{easy_var},sess_stop_window)' < sated_thold & ...
                                  easy_trials > (min_dur * n_trials),1,'First');
    first_unmotivated_trial_t   = session(i).cat_paq(easy_stim).stim_frames{easy_var}(first_unmotivated_trial);
    num_stim_vars = cellfun(@(x) numel(x),{session(i).paq(p).stims(:).vars});
    stims = find(num_stim_vars>0);
    for s = stims
        for v = 1:num_stim_vars(s)
            if ~isempty(first_unmotivated_trial_t)
                flag = session(i).cat_paq(s).stim_frames{v} >= first_unmotivated_trial_t | session(i).cat_paq(s).rxntimes{v} <= min_lick_t;
            else
                flag = session(i).cat_paq(s).rxntimes{v} <= min_lick_t;
            end
            session(i).cat_paq(s).stim_frames{v}(flag) = [];
            session(i).cat_paq(s).responses{v}(flag) = [];
            session(i).cat_paq(s).rxntimes{v}(flag) = [];
        end
    end
    if ~isempty(first_unmotivated_trial_t)
        flag = session(i).stim_var_n_t_rxn_trial(:,4) >= first_unmotivated_trial_t | session(i).stim_var_n_t_rxn_trial(:,5) <= min_lick_t;
    else
        flag = session(i).stim_var_n_t_rxn_trial(:,5) <= min_lick_t;
    end
    session(i).stim_var_n_t_rxn_trial(flag,:) = [];

    % add global stim trial order to cat paq
    for s = 1:numel(session(i).cat_paq)
        num_vars = numel(session(i).cat_paq(s).stim_frames);
        if num_vars>0
            for v = 1:num_vars
                num_trials = numel(session(i).cat_paq(s).stim_frames{v});
                session(i).cat_paq(s).stim_order{v} = zeros(num_trials,1);
                for t = 1:num_trials
                    session(i).cat_paq(s).stim_order{v}(t) = find(session(i).stim_var_n_t_rxn_trial(:,4)==session(i).cat_paq(s).stim_frames{v}(t));
                end
            end
        end
    end
end

% Flag photostimulated cells for each variation
xy_thold        = 10;
plane_depths    = [99 66 33 0];
for i = 1:numel(session)
    % convert z indices into depths for cells and targets
    session(i).s2p.centroids_proc = etlMagTransform(double(session(i).s2p.centroids));
    session(i).overlaps = struct;
    for j = 1:numel(session(i).targets)
        stim = session(i).targets(j).stim_type;
        var = session(i).targets(j).stim_var;
        if ~isempty(session(i).targets(j).targets_xyz)
            session(i).targets(j).targets_xyz_proc = etlMagTransform(session(i).targets(j).targets_xyz);
            % xy distance
            xy_d = euclidean_distance(session(i).s2p.centroids_proc(:,1:2),session(i).targets(j).targets_xyz_proc(:,1:2));
            min_xy = min(xy_d,[],2);
            % flag photostimulated cells
            session(i).overlaps(stim).flag(:,var+1) = min_xy<=xy_thold;
        end
    end  
end

%% Process calcium traces

% trace processing
planes                  = unique([session(1).s2p.stat(:).iplane]);
num_planes              = numel(planes);
p_baseline              = 0.33;
imaging_rate            = session(1).paq(1).frames.rate/num_planes;
imaging_rate_rounded    = round(imaging_rate);
running_baseline_s      = 60;

% main loop
for s = 1:numel(session)
    fprintf(['\nProcessing expt ' num2str(s) ', ' session(s).dir '\n'])
    
    % find stim and bad frames to remove
    session(s).s2p.ops.corrXY_proc      = detrend(session(s).s2p.ops.corrXY) + mean(session(s).s2p.ops.corrXY);
    session(s).s2p.ops.outlierframes    = session(s).s2p.ops.corrXY_proc <= 0;
    session(s).s2p.ops.stimframes       = 0*session(s).s2p.ops.badframes;
    session(s).s2p.ops.stimframes(cell2mat(session(s).cat_paq_raw(7).stim_frames') + round(0:0.6*imaging_rate)) ...
                                        = true;
    
    session(s).s2p.ops.npsc_exclude     = session(s).s2p.ops.stimframes | session(s).s2p.ops.outlierframes | session(s).s2p.ops.badframes;
    session(s).s2p.ops.deconv_exclude   = session(s).s2p.ops.stimframes | session(s).s2p.ops.outlierframes;
    session(s).s2p.ops.artefact_exclude = session(s).s2p.ops.stimframes | session(s).s2p.ops.outlierframes;
    d_bad                               = [0 diff(session(s).s2p.ops.artefact_exclude)];
    session(s).s2p.ops.artefact_starts  = find(d_bad == 1)-1;
    session(s).s2p.ops.artefact_stops   = find(d_bad == -1);
    
    % subtract neuropil, remove artefacts
    session(s).s2p.npsc                 = estimateNeuropilCoefficients(session(s).s2p.f(:,~session(s).s2p.ops.npsc_exclude),session(s).s2p.np(:,~session(s).s2p.ops.npsc_exclude));
    session(s).f_sub                    = halo_subtraction(session(s).s2p.f,session(s).s2p.np,session(s).s2p.npsc,p_baseline);
    session(s).f_sub                    = rm_artifact(session(s).f_sub,session(s).s2p.ops.artefact_starts,session(s).s2p.ops.artefact_stops);

    % calculate whole-trace dff
    session(s).dff                      = bl_normalize_trace(session(s).f_sub,running_baseline_s,imaging_rate_rounded);
    
    % construct run and lick traces
    licks = [];
    running = [];
    for j = 1:numel(session(s).paq)
        licks                           = [licks ; round((session(s).paq(j).licks.frames+1)/num_planes) + session(s).paq(j).cat_offset];
        running                         = [running downsample_HD(session(s).paq(j).running(1:session(s).s2p.ops.frames_per_file(j)*num_planes),num_planes)];
    end
    session(s).run_trace                = running;
    session(s).lick_trace               = 0*session(s).s2p.f(1,:);
    for l = 1:numel(licks)
        session(s).lick_trace(licks(l)) = session(s).lick_trace(licks(l)) + 1;
    end
    session(s).lick_trace_proc          = gaussfilt1d(session(s).lick_trace,imaging_rate_rounded);

end

% remove unused fields
clear f_sub sp licks running
for s = 1:numel(session)
    session(s).s2p = rmfield(session(s).s2p,{'f' 'np' 'spks'});
end


%% Return trial-wise data in [sessions trial_type] cell array
% stas
pre_frames_sta.ps       = round(1*imaging_rate);
post_frames_sta.ps      = round(4*imaging_rate);
stim_frame.ps           = pre_frames_sta.ps+1;
pre_window_dff.ps       = 1:pre_frames_sta.ps;
post_window_dff.ps      = stim_frame.ps + ceil((0.6*imaging_rate):(1*imaging_rate));
stims                   = [6 7];
nogo_stim               = 6;
go_stim                 = 7;

[stim_vars,num_cells] = deal([]);
targets_sv = [[session(1).targets(:).stim_type] ; [session(1).targets(:).stim_var]+1];
for stim = 1:numel(session(1).cat_paq)
    if ~isempty(session(1).cat_paq(stim).stim_frames)
        for var = 1:numel(session(1).cat_paq(stim).stim_frames)
            stim_vars = [stim_vars [stim ; var]];
            num_cells = [num_cells size(session(1).targets(min(targets_sv == stim_vars(:,end),[],1)).targets_xyz,1)];
        end
    end
end
num_stim_vars = size(stim_vars,2);
num_sessions = numel(session);
thold_var               = num_cells == 50;
ng_var                  = num_cells == 0;
max_var                 = num_cells == 200;

[tw_traces_raw,tw_traces_dff,mean_traces_dff,tw_response_dff,...
    mean_response_dff,tw_traces_lick,mean_traces_lick,tw_traces_run,...
    mean_traces_run,tw_response,tw_rxn,...
    tw_traces_zs,mean_traces_zs,tw_response_zs,mean_response_zs,...
    targeted_flag] ...
    = deal(cell(num_sessions,num_stim_vars));
for sess = 1:num_sessions
    fprintf(['\nCalculating trial data for expt ' num2str(sess) ', ' session(sess).dir '\n'])
    
    for trial_type = 1:num_stim_vars
        % make stas [neurons trials frames]
        stim                                        = stim_vars(1,trial_type);
        var                                         = stim_vars(2,trial_type);
        these_frames                                = session(sess).cat_paq(stim).stim_frames{var};
        
        % find targets on each trial
        if ~isempty(session(sess).overlaps(stim).flag)
            targeted_flag{sess,trial_type}          = session(sess).overlaps(stim).flag(:,var);
        else
            targeted_flag{sess,trial_type}          = false(size(session(sess).f_sub,1),1);
        end
        
        % raw f, dff
        tw_traces_raw{sess,trial_type}              = make_sta_traces(session(sess).f_sub,these_frames,pre_frames_sta.ps,post_frames_sta.ps);
        bl                                          = nanmean(tw_traces_raw{sess,trial_type}(:,:,pre_window_dff.ps),3);
        tw_traces_dff{sess,trial_type}              = (tw_traces_raw{sess,trial_type} - bl) ./ bl;
        mean_traces_dff{sess,trial_type}            = squeeze(nanmean(tw_traces_dff{sess,trial_type},2));
        tw_response_dff{sess,trial_type}            = mean(tw_traces_dff{sess,trial_type}(:,:,post_window_dff.ps),3);
        mean_response_dff{sess,trial_type}          = mean(tw_response_dff{sess,trial_type},2);
        
        % z-scored response
        bl_sd                                       = nanstd(tw_traces_raw{sess,trial_type}(:,:,pre_window_dff.ps),[],3);
        bl_mu                                       = nanmean(tw_traces_raw{sess,trial_type}(:,:,pre_window_dff.ps),3);
        tw_traces_zs{sess,trial_type}               = (tw_traces_raw{sess,trial_type} - bl_mu) ./ bl_sd;
        mean_traces_zs{sess,trial_type}             = squeeze(mean(tw_traces_zs{sess,trial_type},2));
        tw_response_zs{sess,trial_type}             = mean(tw_traces_zs{sess,trial_type}(:,:,post_window_dff.ps),3);
        mean_response_zs{sess,trial_type}           = mean(tw_response_zs{sess,trial_type},2);
        
        % licks
        tw_traces_lick{sess,trial_type}             = squeeze(make_sta_traces(session(sess).lick_trace,these_frames,pre_frames_sta.ps,post_frames_sta.ps));
        mean_traces_lick{sess,trial_type}           = nanmean(tw_traces_lick{sess,trial_type},1);
        
        % running
        tw_traces_run{sess,trial_type}              = squeeze(make_sta_traces(session(sess).run_trace,these_frames,pre_frames_sta.ps,post_frames_sta.ps));
        mean_traces_run{sess,trial_type}            = nanmean(tw_traces_run{sess,trial_type},1);
        
        % trial info
        tw_response{sess,trial_type}                = session(sess).cat_paq(stim).responses{var}==1;
        tw_rxn{sess,trial_type}                     = session(sess).cat_paq(stim).rxntimes{var};
    end
end

%% Spontaneous lick analysis
resp_threshold      = 1;                                % z-scored
min_bout_len        = 3;                                % minimum no. licks in a bout
min_d               = round(imaging_rate_rounded/2);    % minimum interval between consecutive licks for them to be included in the same bout
pre_buffer          = imaging_rate_rounded;             % only look for spontaneous licks separated from preceding licks by this amount
trial_epoch         = 0:round(4*imaging_rate);          % trial-periods to exclude from spontaneous licking
maxLag              = imaging_rate_rounded*8;
cc_len              = maxLag*2+1;

[bout_starts,bout_n] = deal(cell(num_sessions,1));
[tw_lick_traces_dff,mean_lick_traces_dff,tw_lick_response_dff,...
    mean_lick_response_dff,lick_response_reliability,lick_response_cc,...
    tw_lick_traces_zs,mean_lick_traces_zs,tw_lick_response_zs,mean_lick_response_zs,...
    lick_exclude,stim_trace,spont_lick_trace,spont_lick_cc,spont_lick_p,...
    tw_spont_lick_traces,mean_spont_lick_traces,...
    crossCorrDff,crossCorrSign,crossCorrDffPosMean,crossCorrDffNegMean,...
    staTracesDffPosMean,staTracesDffNegMean] ...
    = deal(cell(num_sessions,1));

[lick_cc_tholds] = deal(zeros(num_sessions,1));
for s = 1:num_sessions
    fprintf(['\nAnalysing lick responses for expt ' num2str(s) ', ' session(s).dir '\n'])
    
    % Return stim and lick traces
    stim_frames                 = session(s).stim_var_n_t_rxn_trial_raw(session(s).stim_var_n_t_rxn_trial_raw(:,1)==7,4);
    trial_periods               = stim_frames + trial_epoch;
    stim_trace{s}               = 0*session(s).lick_trace;
    stim_trace{s}(trial_periods)= true;
    lick_trace                  = session(s).lick_trace;
    
    % Find lick bouts
    [bout_starts{s},bout_n{s}]  = find_lick_bouts(find(lick_trace),min_d,pre_buffer);
    bout_stim                   = stim_trace{s}(bout_starts{s});
    flag                        = bout_n{s}>=min_bout_len & ~bout_stim;
    bout_starts{s}              = bout_starts{s}(flag);
    bout_n{s}                   = bout_n{s}(flag);
    
    % Trace correlation with licking
    spont_lick_trace{s}         = double(gaussfilt1d(session(s).lick_trace(~stim_trace{s}),round(imaging_rate_rounded*0.5)));
    dff                         = double(session(s).dff(:,~stim_trace{s}));
    lt                          = find(spont_lick_trace{s});
    lick_epoch                  = lt(1):lt(end);
    spont_lick_trace{s}         = spont_lick_trace{s}(lick_epoch);
    dff                         = dff(:,lick_epoch);
    [spont_lick_cc{s},spont_lick_p{s}] ...
                                = corr(dff',spont_lick_trace{s}');
    
    % post lick STA (dff)
    [~,tw_lick_traces_dff{s},mean_lick_traces_dff{s}] ...
                                = make_sta_traces_dff(session(s).f_sub,bout_starts{s},pre_frames_sta.ps,post_frames_sta.ps);
    tw_lick_response_dff{s}     = mean(tw_lick_traces_dff{s}(:,:,pre_window_dff.ps(end)+[1:imaging_rate_rounded]),3);
    mean_lick_response_dff{s}   = nanmean(tw_lick_response_dff{s},2);
    for i = 1:size(tw_lick_traces_dff{s},1)
        lick_response_cc{s}(i,1) ...
                                = nanmean(nanmean(off_diag(corr(squeeze(tw_lick_traces_dff{s}(i,:,:))'))));
    end
    
    % post lick STA (zs)
    [~,tw_lick_traces_zs{s},mean_lick_traces_zs{s}] ...
                                = make_sta_traces_zs(session(s).f_sub,bout_starts{s},pre_frames_sta.ps,post_frames_sta.ps);
    tw_lick_response_zs{s}      = mean(tw_lick_traces_zs{s}(:,:,pre_window_dff.ps(end)+[1:imaging_rate_rounded]),3);
    mean_lick_response_zs{s}    = nanmean(tw_lick_response_zs{s},2);
    lick_response_reliability{s} = max([mean(tw_lick_response_zs{s}>resp_threshold,2) mean(tw_lick_response_zs{s}<-resp_threshold,2)],[],2);
    
    tmp                         = lick_response_cc{s};
    lick_cc_tholds(s)           = abs(prctile(tmp(tmp<0),5));
    
    % Exclude licking cells
    lick_exclude{s}             = spont_lick_p{s} < 0.05 | ...
                                  abs(mean_lick_response_zs{s}) > resp_threshold & lick_response_reliability{s} > 0.25;%| lick_response_cc{s} > lick_cc_tholds(s);
    
    % lick STA (licking trace)
    tw_spont_lick_traces{s}     = make_sta_traces_dff(session(s).lick_trace,bout_starts{s},pre_frames_sta.ps,post_frames_sta.ps);
    tw_spont_lick_traces{s}     = squeeze(tw_spont_lick_traces{s});
    mean_spont_lick_traces{s}   = nanmean(tw_spont_lick_traces{s},1);
                                       
    % Cross-correlation
    crossCorrDff{s}             = zeros(cc_len,size(dff,1));
    zLick                       = zscore(spont_lick_trace{s});
    zdff                        = zscore(dff,[],2);
    for i = 1:size(dff,1)
        [crossCorrDff{s}(:,i),lags]...
                                = xcorr(zdff(i,:)',zLick',maxLag,'coeff');
    end
    crossCorrDff{s}             = crossCorrDff{s}';
    [~,m]                       = max(abs(crossCorrDff{s}),[],2);
    idx                         = sub2ind(size(crossCorrDff{s}),1:size(crossCorrDff{s},1),m');
    crossCorrSign{s}            = sign(crossCorrDff{s}(idx))';
    crossCorrDffPosMean{s}      = crossCorrDff{s}(crossCorrSign{s}>0 & lick_exclude{s},:);
    crossCorrDffNegMean{s}      = crossCorrDff{s}(crossCorrSign{s}<0 & lick_exclude{s},:);
    staTracesDffPosMean{s}      = mean_lick_traces_dff{s}(crossCorrSign{s}>0 & lick_exclude{s},:);
    staTracesDffNegMean{s}      = mean_lick_traces_dff{s}(crossCorrSign{s}<0 & lick_exclude{s},:);
end

