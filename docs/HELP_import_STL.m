%% import_STL
% Below is a demonstration of the features of the |import_STL| function

%% Syntax
% |[stlStruct] = import_STL(fileName);|

%% Description
% Use |import_STL| to import binary or txt type STL files.

%% Examples

clear; close all; clc;

%%
% Plot settings
fontSize=25; 

%%
% Set path and file name

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 

%File name for STL
fileName=fullfile(pathName,'vertebra.stl'); 

%% Import a txt type STL file as patch data
[stlStruct] = import_STL(fileName);

%%
% Access the data from the STL struct
F=stlStruct.solidFaces{1}; %Faces
V=stlStruct.solidVertices{1}; %Vertices

%%
% In STL files nodes are not shared between triangles. Therefore the need
% to be merged with they are intended to be shared. 

[F,V]=mergeVertices(F,V); % Merging nodes

%%
% Plotting the model

cFigure;
title('Imported patch data from STL','fontSize',25);
gpatch(F,V,'gw');
axisGeom;
camlight('headlight');
lighting phong; axis off; 
gdrawnow;

%% Import a binary type STL file as patch data

fileName=fullfile(pathName,'femur_bin.stl'); 
[stlStruct] = import_STL(fileName);

%%
% Access the data from the STL struct
F=stlStruct.solidFaces{1}; %Faces
V=stlStruct.solidVertices{1}; %Vertices

%%
% Merge nodes
[F,V]=mergeVertices(F,V); % Merging nodes

%%
% Plotting the model

cFigure;
title('Imported patch data from STL','fontSize',fontSize);
gpatch(F,V,'gw');
axisGeom;
camlight('headlight');
lighting phong; axis off; 
gdrawnow;

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
