function [bhv] = pyb_get_stim_powers(bhv,path)
% finds LED power files in path for the trials in the session passed via
% bhv

% load all var files for this date
thisDate = strsplit(bhv(1).trialstart,'_');
[varFiles,varFileNames] = pyb_load_varfiles(path,thisDate{1});

% if there are varFiles
if ~isempty(varFiles)
    % for each trial
    for t = 1:numel(bhv)
        % only stim trials
        if bhv(t).stim_type == 5
            
            if ~isempty(bhv(t).trialstart)
                % find trial start
                trialstartString = bhv(t).trialstart;
                % find relevant var file for this trial
                thisVarFileName = pyb_return_varfile(path,trialstartString);
            else
                thisVarFileName = [];
            end
            
            % if one exists and the number of variations matches,
            % find the relevant power for this variation else set
            % power to nan
            numVars = numel(unique([bhv(:).stim_var]));
            if ~isempty(thisVarFileName) ...
                    && numVars <= numel(varFiles{strcmp(varFileNames,thisVarFileName)}.pulse_powers)
                thisVar = bhv(t).stim_var+1;
                thisPower = varFiles{strcmp(varFileNames,thisVarFileName)}.pulse_powers(thisVar);
                bhv(t).power = thisPower;
            else
                bhv(t).power = nan;
            end
        else
            bhv(t).power = 0;
        end
    end
    % if no varFiles, set all stim trial powers to nan
else
    for t = 1:numel(bhv)
        if bhv(t).stim_type == 5
            bhv(t).power = nan;
        else
            bhv(t).power = 0;
        end
    end
end

end

