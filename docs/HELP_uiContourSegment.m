%% uiContourSegment
% Below is a demonstration of the features of the |uiContourSegment| function

%%
clear; close all; clc;

%% Loading example image data
% Loading MRI data
toolboxPath=fileparts(fileparts(mfilename('fullpath')));

imageFolder=fullfile(toolboxPath,'data','DICOM','KNEE_UTE','IMDAT');

loadName=fullfile(imageFolder,'IMDAT.mat');
IMDAT_set=load(loadName);
M=double(IMDAT_set.type_1);
M_info=IMDAT_set.type_1_info(1);
[v,OR,r,c]=dicom3Dpar(M_info);

%% Visualizing data

hf=sliceViewer(M,v,2);

%% REMOVING BACKGROUND THROUGH THRESHOLDING/DILATION PROCEDURE
% Removal of background may be useful see |uiThreshErode|
L_BG=true(size(M)); %In this example background is not removed

thresholdInitial=0.1; %with respect to normalised image
preBlurKernalSize=0; %with respect to normalised image
groupCropOption=0;
% [L_BG]=uiThreshErode(M,thresholdInitial,preBlurKernalSize,groupCropOption);

%% SETTING CONTROL PARAMETERS
cPar.minContourSize=150;            %Minimal size of detected contour
cPar.smoothFactor=0.15;              %Degree of smoothing csaps function (cubic smoothing spline)
cPar.pointReductionFactor=1;        %Reduction factor for contour smoothening
cPar.logicBackGround=L_BG;   %Ones (white) describe image data regions of interest i.e. a mask
cPar.v=v;                           %Voxel size
cPar.recoverOn=0; %Turn on or off file recovery mode 
cPar.sliceRange=98:144; %This can be a custom range. For unvisited slices the contour is empty
saveName=[];%'fibula';                        %If not empty this is where the contours are saved

%% SEGMENTING CONTOURS
% Run the following code:
% |[Vcs]=uiContourSegment(M,cPar,saveName)|;
% Use the various options to try to segment the outer skin contour, press H
% to reveal the help window. 

%%
% Uncomment to test
[Vcs]=uiContourSegment(M,cPar,saveName);

%% 
% Example contour process: 
% 
% <<contour_1.png>>

%%
% Clearly only features "detectable" with contour levels can be segmented or
% otherwise manual drawing is needed.

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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
