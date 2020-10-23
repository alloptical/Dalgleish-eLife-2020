function [sessIdx] = pyb_file2id(bhvFileStr)
% walks through cell array of bhv file names (animal_date_time) to and
% returns an array indexing the unique dates (i.e. all files from a given
% date will have the same number)

currDate = 0;
sessCount = 0;
sessIdx = zeros(numel(bhvFileStr),1);
for i = 1:numel(bhvFileStr)
    [~,fl] = fileparts(bhvFileStr{i});
    tmp = strsplit(fl,'_');
    thisDate = str2num(tmp{2});
    if thisDate ~= currDate
        sessCount = sessCount+1;
        currDate = thisDate;
    end
    sessIdx(i) = sessCount;
end

end

