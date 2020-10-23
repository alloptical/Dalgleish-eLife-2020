function [catBhvData,catStimVars,catFilenames] = pyb_cat_daily_sessions(bhvData,fileNames)
% combines bhv sessions from the same day

sessIdx = pyb_file2id(fileNames);
uniqueIdcs = unique(sessIdx);
nSess = numel(uniqueIdcs);
[catBhvData,catStimVars,catFilenames] = deal(cell(nSess,1));
for i = 1:numel(uniqueIdcs)
    theseIdcs = find(sessIdx==uniqueIdcs(i));
    for j = 1:numel(theseIdcs)
        catBhvData{i} = [catBhvData{i} ; bhvData{theseIdcs(j)}];
        catFilenames{i} = [catFilenames{i} ; fileNames(theseIdcs(j))];
        sv = pyb_return_stimvars(bhvData{theseIdcs(j)});
        if isempty(catStimVars{i})
            catStimVars{i} = sv;
        else 
            toAdd = false(1,size(sv,2));
            for k = 1:size(sv,2)
                toAdd(k) = all(any(catStimVars{i}~=sv(:,k),1));
            end
            catStimVars{i} = [catStimVars{i} sv(:,toAdd)];
        end
    end
end

end

