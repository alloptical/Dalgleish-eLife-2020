function [] = unzipDirectories(zipDir)
%% Unzip directories
% This function will walk through a directory of zipped directorys and
% unzip them all. As it goes it will check to make sure that the
% original directory structure of each .zip file is preserved.

% This sounds trivial, but sometimes Windows will add an additional
% directory, of the same name, upstream of the zipped directory, e.g.:
% myDir/ --> myDir.zip --> myDir/myDir/
%
% Using this function ensures that: 
% myDir/ --> myDir.zip --> myDir/
%
% Input: directory containing .zip files of directories you want to unzip
%
% Henry Dalgleish 2020

%%
d = dir(zipDir);
zipFiles = d(cellfun(@(x) contains(x,'.zip'),{d(:).name}));
if numel(zipFiles)>0
    cd(zipDir);
    for i = 1:numel(zipFiles)
        pathfile = fullfile(zipFiles(i).folder,zipFiles(i).name);
        % unzip files
        unzip(pathfile);
        % if on Windows need to account for the fact that dir.zip will get
        % unzipped to dir/dir/ instead of dir/ (like on mac)...
        [p,n,e] = fileparts(pathfile);
        if exist(fullfile(p,n,n),'dir') == 7
            movefile(fullfile(p,n),fullfile(p,'tmp')) % rename upstream dir
            movefile(fullfile(p,'tmp',n),fullfile(p,n)) % move downstream dir upstream
            rehash()
            rmdir(fullfile(p,'tmp'),'s') % remove renamed upstream dir
        end
    end
end

end

