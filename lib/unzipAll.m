function unzipAll(varargin)

% function unzipAll(pathName,pathNameOutput)
% ------------------------------------------------------------------------
%
% This function unzips all zip files in the folder pathName. The first
% input is the folder path name, the following inputs are possible inputs
% of the unzip function. 
%
% See also: dicomanon
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/09: Created
%------------------------------------------------------------------------

%%
switch nargin 
    case 1
        pathName=varargin{1};
        varOutput=pathName;
    case 2
        pathName=varargin{1};
        varOutput=varargin{2};
end
%%

if ischar(varOutput) %Process only input folder
    files = dir(fullfile(pathName,'*.zip'));
    files={files(1:end).name};
    files=sort(files(:));
    if ~isempty(files)
        if numel(files)>1
            hw = waitbar(0,'Unzipping files...');
            for q=1:1:numel(files)
                fileName=fullfile(pathName,files{q});
                unzip(fileName,varOutput);
                waitbar(q/numel(files),hw,['Unzipping files...',num2str(round(100.*q/numel(files))),'%']);
            end
            close(hw)
        else
            fileName=fullfile(pathName,files{1});
            unzip(fileName,varOutput);
        end
    end
else %Process input folder and all subfolders
    %Process current
    unzipAll(pathName,pathName); 
    
    %Process sub-folders
    [pathNames]=getSubPaths(pathName);
    if ~isempty(pathNames)
        for q=1:1:numel(pathNames)
            pathNameNow=pathNames{q};
            unzipAll(pathNameNow,pathNameNow);
        end
    end
end






