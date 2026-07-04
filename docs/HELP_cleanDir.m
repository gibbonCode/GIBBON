%% cleanDir
% Below is a demonstration of the features of the |cleanDir| function

%%
clear; close all; clc;

%% Syntax
% |cleanDir(pathName,extCell);|

%% Description 
% This function can clean a folder by removing either all files in a folder
% or all files with the extensions specified in the cell array extCell.
%
% E.g. if one is interested in removing all .txt and .xml files from the
% folder specified by pathName on could use: 
% pathName='D:\MATLAB\GIBBON'
% extCell={'txt','xml'}; %Extensions of files to delete
% cleanDir(pathName,extCell);

%% Examples 
% 

%%
% Create an example folder containing files with 3 different extensions,
% i.e. txt, csv and doc. 

% Create temporary example folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','temp','cleanCheck');
mkdir(pathName)

% Create files
n=3; %Number of files for each type
for q=1:1:n %Create n files
    fileID=fopen(fullfile(pathName,['temp',num2str(q),'.txt']),'w');
    fprintf(fileID,'%d\n',pi);
    fclose(fileID);
end

for q=1:1:n %Create n files
    fileID=fopen(fullfile(pathName,['temp',num2str(q),'.csv']),'w');
    fprintf(fileID,'%d\n',pi);
    fclose(fileID);
end

for q=1:1:n %Create n files
    fileID=fopen(fullfile(pathName,['temp',num2str(q),'.xml']),'w');
    fprintf(fileID,'%d\n',pi);
    fclose(fileID);
end

%%
% Show current folder content
disp('Old folder content:')
ls(pathName)

%%
% Clean up directory by removing all files whose extensions are a member of
% the extension set provided. 

extCell={'csv','txt'}; %File extensions for files to remove
cleanDir(pathName,extCell)

%%
% Show current folder content
disp('New folder content:')
ls(pathName)

%%
% remove the temporary folder created for this example
rmdir(pathName,'s')

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
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
