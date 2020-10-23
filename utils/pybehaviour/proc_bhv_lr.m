function [bhv,params] = proc_bhv_lr(bhv,varargin)
% Lower level behaviour processing function

params.sated_thold = 0.85;
params.min_dur = 0.5;
params.min_lick_t = 0.3;
params.response_window = 1;
params.remove_ambiguous = true;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'sated_thold')
        params.sated_thold = varargin{v+1};
    elseif strcmpi(varargin{v},'min_dur')
        params.min_dur = varargin{v+1};
    elseif strcmpi(varargin{v},'min_lick_t')
        params.min_lick_t = varargin{v+1};
    elseif strcmpi(varargin{v},'response_window')
        params.response_window = varargin{v+1};
    elseif strcmpi(varargin{v},'remove_ambiguous')
        params.remove_ambiguous = varargin{v+1};
    end
end

% find final good trial
easy_trials = find([bhv(:).stim_type]==7 & [bhv(:).stim_var]==0);
easy_responses = movmean([bhv(easy_trials).responded],10);
sated_trials = easy_trials(easy_responses<params.sated_thold & easy_trials>(numel(bhv)*params.min_dur));
if isempty(sated_trials)
    final_good_trial = numel(bhv);
else
    final_good_trial = sated_trials(1)-1;
end

% remove sated trials, early trials
bhv = bhv(1:final_good_trial);
to_keep = [bhv(:).rxntime]>params.min_lick_t | isnan([bhv(:).rxntime]);
bhv = bhv(to_keep);

% correct response window, remove autorewarded trials
toremove = false(numel(bhv),1);
for i = 1:numel(bhv)
    bhv(i).parameters.responseWindow = params.response_window;
    bhv(i).parameters.params = params;
    if ~isnan(bhv(i).rxntime)
        % trials with responses > responsWindow are misses
        if bhv(i).rxntime > bhv(i).parameters.responseWindow
            bhv(i).rxntime = nan;
            bhv(i).responded = false;
        % auto-rewarded trials with responses in responseWindow but
        % following auto-reward time by sufficiently long time interval are
        % removed (ambiguous)
        elseif bhv(i).rxntime <= bhv(i).parameters.responseWindow ...
                && bhv(i).auto_reward == 1 ...
                && bhv(i).rxntime >= (bhv(i).parameters.autoRewardDelay + params.min_lick_t) ...
                && bhv(i).response_required ~= 0
            if params.remove_ambiguous
                toremove(i) = true;
            else
                bhv(i).responded = false;
            end
        end
    end
end
bhv(toremove) = [];

end

