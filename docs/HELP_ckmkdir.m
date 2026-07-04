%% ckmkdir
% Below is a demonstration of the features of the |ckmkdir| function

%%
clear; close all; clc;

%% Syntax
% |[existFlag]=ckmkdir(dirName);|

%% Description 
% Check for the existance of the directory dirName and make the directory
% if it does not. The optional output is the exitFlag which is 0 if it did
% not exist and 1 if it did. 
%
% See also |mkdir|

%% Examples 
% 

%% Make a directory which does not exist yet

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','temp');

% % Make new directory first
% mkdir(fullfile(pathName,'newDir'))

disp('Old folder content: ')
ls(pathName)

%%
% Create directory if not present

% File name 
dirName=fullfile(pathName,'new_folder');

[existFlag]=ckmkdir(dirName)

disp('New folder content: ')
ls(pathName)

%%
% If the directory is already there the folder is ignored and existFlag is
% true 

[existFlag]=ckmkdir(dirName)

%%
% remove the temporary folder created for this example
rmdir(dirName,'s')


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
