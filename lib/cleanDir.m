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
    
    for q=1:1:numel(fileList) %For all files
        fileName=fullfile(pathName,fileList{q}); %Current file name
        if ~isfolder(fileName) %If it is not a directory
            delete(fileName); %Delete the file
        end
    end
else %Delete files matching extension list
    for qc=1:1:numel(extCell) %For all extensions
        ext=extCell{qc}; %Current extension
        fileList = dir(fullfile(pathName,['*.',ext]));
        fileList={fileList(1:end).name}; %Current file list
        
        if ~isempty(fileList) %For all files matching the extenion            
            for q=1:1:numel(fileList)
                fileName=fullfile(pathName,fileList{q}); %Current file name
                delete(fileName); %Delete the file
            end
        end
    end
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
