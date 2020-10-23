function [thisVar] = pyb_return_varfile(path,trialstartString)

% find varFiles generated on this date
trialDateTime = strsplit(trialstartString,'_');
varFiles = return_fullfile(path,'VarFile',trialDateTime{1});
if ~iscell(varFiles)
    varFiles = {varFiles};
end

% find all varFile timestamps and sort in ascending order
varTimes = cellfun(@(x) str2num(myStrsplit(myFileparts(x,'name'),'_',4)),varFiles);
[varTimes,order] = sort(varTimes,'Ascend');
varFiles = varFiles(order);
trialTime = str2num(trialDateTime{2});

% find most recent varFile relative to this trial
thisIdx = find(varTimes<=trialTime,1,'Last');
if ~isempty(thisIdx)
    thisVar = varFiles{thisIdx};
else
    thisVar = [];
end

end

