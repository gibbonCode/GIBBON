%% plotContours
% Below is a demonstration of the features of the |plotContours| function

%%
clear; close all; clc;

%% Syntax
% |[handleCell]=plotContours(contourSet,optionStruct);|

%% Description 
% This function plots segmented contours 

%% Import image data for this demo
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
pathNameImageData=fullfile(defaultFolder,'data','DICOM','0001_human_calf');
loadNameImageData=fullfile(pathNameImageData,'IMDAT','IMDAT.mat');
IMDAT_struct=load(loadNameImageData); %The image data structure
G = IMDAT_struct.G; %Geometric/spatial information
v=G.v; %The voxel size
M= IMDAT_struct.type_1; %The image data

%% Example: Visualize a segmented contour  
% In this example the contour data is first loaded and then the plotting is
% done. 

%%
% Visualize image data
sv3(M,v);
drawnow;

%%
% Add contour plots 
contourName='imseg_calf_tibia'; %Contour name
optionStruct.pathName=fullfile(defaultFolder,'data','imseg'); %Folder name for contours

%Plot contours
plotContours(contourName,optionStruct);

%% Example: Visualize a set of segmented contours  
% In this example the contour data is loaded by the plotting function. Only
% contour file and path names are passed to the plot function. 

%%
% Visualize image data
sv3(M,v);
drawnow;

%%
% Add contour plots 

%Contour names
contourSet={'imseg_calf_skin',...
            'imseg_calf_muscle_fat',...
            'imseg_calf_tibia',...
            'imseg_calf_fibula'}; 

%Fill option structure
optionStruct.LineWidth=5; %Line width
optionStruct.Color=gjet(numel(contourSet));     
optionStruct.pathName=fullfile(defaultFolder,'data','imseg'); %Folder name for contours
optionStruct.hAxis=gca; %Axis to plot in

%Plot contours
handleCell=plotContours(contourSet,optionStruct);

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
