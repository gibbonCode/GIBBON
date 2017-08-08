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
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
