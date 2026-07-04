%% getInnerVoxel
% Below is a demonstration of the features of the |getInnerVoxel| function

%%
clear; close all; clc;

%% Syntax
% |[indInternal]=getInnerVoxel(L,searchRadius,plotOn);|

%% Description 
% Obtains a voxel index for a voxel that is part of the interior of the
% labelled image L. searchRadius defines the approximate distance to the
% boundary. plotOn visualizes the label image and interior point.

%% Examples 
% 
%%
% Plot settings
fontSize=10;
faceAlpha1=1;
faceAlpha2=0.3;
cMap=[0.5 0.5 0.5; gjet(4)];

%%
% 
% Create example image
[F,V]=graphicsModels(8);

[M,G,~]=patch2Im(F,V,[],0.05);
L=~isnan(M);

%%
% Get index for an interior voxel
searchRadius=3;
plotOn=1;
[indInternal]=getInnerVoxel(L,searchRadius,plotOn)

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
