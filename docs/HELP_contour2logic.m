%% contour2logic
% Below is a demonstration of the features of the |contour2logic| function

%%
clear; close all; clc;

%% Syntax
% |[varargout]=contour2logic(M,v,Vcs);|

%% Description 
% This function converts contours to a logic or labelled data. The logics
% represent wether voxels are in, on, or outside the contour. 
%
% The input consists of: 
%%
% 
% * A 3D image |M| (or alternatively the size of M). 
% * The |vozelSize| a 1x3 vector specifying the size of the voxels in the
% row, column, and slice direction.
% * A cell array |Vcs| containing one or more contours per slice. If the
% image has n slices then Vcs should be an nx1 cell array, i.e. contours
% are defined for each slice. 

%% Examples 
% 
%% Import image data for this demo

defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
pathNameImageData=fullfile(defaultFolder,'data','DICOM','0001_human_calf');
loadNameImageData=fullfile(pathNameImageData,'IMDAT','IMDAT.mat');
IMDAT_struct=load(loadNameImageData); %The image data structure
G = IMDAT_struct.G; %Geometric/spatial information
v=G.v; %The voxel size
M= IMDAT_struct.type_1; %The image data

%%
contourName='imseg_calf_tibia';
pathName=fullfile(defaultFolder,'data','imseg'); %Folder name for contours

%% Compute levelset

loadName=fullfile(pathName,contourName);
load(loadName); %Load segmentation structure
Vcs=saveStruct.ContourSet; %Access the contour data

[logicIn,logicOn,N]=contour2logic(M,v,Vcs);

%%
% Visualize logic image and contours together

%Visualize logic image
sv3(logicIn,v); %Open slice viewer for levelset

%Visualize contours
optionStruct.Color='r';
plotContours({Vcs},optionStruct);  %Plot contours

%Add colorbar labels
[~,hc]=icolorbar;
hc.TickLabels={'Out','In'};

drawnow;

%%
% Visualize labelled image and contours together

%Visualize label image
vizStruct.colormap=viridis(3); %Set colormap for levelset visualization
hf2=sv3(N,v,vizStruct); %Open slice viewer for levelset

%Visualize contours
optionStruct.Color='r';
plotContours({Vcs},optionStruct);  %Plot contours

%Add colorbar labels
[~,hc]=icolorbar;
hc.TickLabels={'Out','On','In'};

drawnow;

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
