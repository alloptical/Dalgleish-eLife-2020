function [varFiles,varFileNames] = pyb_load_varfiles(path,date)
% return all LED power files in path for this date

if isempty(date)
    varFileNames = return_fullfile(path,'VarFile');
else
    varFileNames = return_fullfile(path,'VarFile',date);
end
if ~isempty(varFileNames)
    if ~iscell(varFileNames)
        varFileNames = {varFileNames};
    end
    varFiles = cell(numel(varFileNames),1);
    for v = 1:numel(varFileNames)
        varFiles{v} = load(varFileNames{v});
        varFiles{v} = varFiles{v}.o;
    end
else
    varFiles = {};
    varFileNames = {};
end

end

