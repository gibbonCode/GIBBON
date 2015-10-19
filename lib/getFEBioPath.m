function FEBioPath=getFEBioPath

% function FEBioPath=getFEBioPath
% ------------------------------------------------------------------------
% This function reads a user specified FEBio path name from a config text
% file. The text file may contain multiple paths (e.g. for multiple
% operational systems), however the first valid patch is always used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/10/05 Updated for GIBBON
% 2015/10/05 Enabled searching valid paths in multi-line config. file.
%------------------------------------------------------------------------


filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths
try 
    [T]=txtfile2cell(fileName);    
catch
    T={};
end

%Get first valid path
FEBioPath=[];
if ~isempty(T)
    for q=1:1:numel(T)
        pathName=T{q};
        if exist(pathName,'file')==2
            FEBioPath=pathName;
            return
        end
    end
end

