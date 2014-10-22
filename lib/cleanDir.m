function cleanDir(varargin)

% function cleanDir(pathName,extCell)
% ------------------------------------------------------------------------
% This function can clean a folder by removing either all files in a folder
% or all files with the extensions specified in the cell array extCell.
%
% E.g. if one is interested in removing all .txt and .xml files from the
% folder specified by pathName on could use: 
% pathName='D:\MATLAB\GIBBON'
% extCell={'txt','xml'}; %Extensions of files to delete
% cleanDir(pathName,extCell);
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------
     
%% Parse input
switch nargin
    case 1
        pathName=varargin{1};
        extCell={};
    case 2
        pathName=varargin{1};
        extCell=varargin{2};
end

%% Removing files from folder

if isempty(extCell) %Delete all files
    fileList = dir(pathName);
    fileList={fileList(1:end).name}; %Current file list
    
    for q=1:1:numel(fileList); %For all files
        fileName=fullfile(pathName,fileList{q}); %Current file name
        if ~isdir(fileName); %If it is not a directory
            delete(fileName); %Delete the file
        end
    end
else %Delete files matching extension list
    for qc=1:1:numel(extCell) %For all extensions
        ext=extCell{qc}; %Current extension
        fileList = dir(fullfile(pathName,['*.',ext]));
        fileList={fileList(1:end).name}; %Current file list
        
        if ~isempty(fileList); %For all files matching the extenion            
            for q=1:1:numel(fileList);
                fileName=fullfile(pathName,fileList{q}); %Current file name
                delete(fileName); %Delete the file
            end
        end
    end
end