%% getInnerPoint
% Below is a demonstration of the features of the |getInnerPoint| function

%%
clear; close all; clc;

%% Syntax
%|[V_inner,M,G,ML,voxelSize]=getInnerPoint(F,V,searchRadius,voxelSize,plotOn)|
% |[varargout]=getInnerPoint(varargin);|

%% Description 
% This function computes an arbitrary interior point for the input
% geometry. The function uses |patch2Im| to convert the geometry to an
% image description. Next the interior voxel set is convolution with a
% spherical kernel of a desired size. An interior point is then chosen
% based on the location with the maximum output in the convoluted image. In
% plainer English this means that an attempt is made to find a point that
% is inside the geometry and approximately the spherical kernel radius
% offset inwards from the boundary. 
% Input consists of the faces F, the vertices V, the searchRadius (kernel
% radius), the voxelSize, and a plotting option plotOn.
% The voxel size should be shall enough such that interior (and not just
% boundary) voxels can be found. Interior voxels and fully inside the
% geometry and do not touch the boundary.  

%%
% Plot settings
markerSize=50; 

%% Examples 

%% Example 1: Basic use to find an arbitrary point inside the input geometry

%%
% Create example geometry
testCase=3;
switch testCase
    case 1
        [F,V]=geoSphere(1,1);
    case 2
        [F,V]=stanford_bunny;
    case 3
        [F,V]=graphicsModels(4);
end

%%
% Find interior point using default settings
V_in=getInnerPoint(F,V);

%%

cFigure; 
subplot(1,2,1); hold on;
gpatch(F,V,'gw');
axisGeom; camlight headlight;

subplot(1,2,2); hold on;
gpatch(F,V,'gw','none',0.5);
plotV(V_in,'r.','MarkerSize',markerSize)
axisGeom; camlight headlight;
drawnow;

%% Example: using full input/output set

D=patchEdgeLengths(F,V); %Get edge lengths
voxelSize=mean(D)/2; %Set voxel size as half of the mean edge length
searchRadius=3*voxelSize; %Use 3 voxel search radius
plotOn=1;
[V_in,M,G,ML,voxelSize]=getInnerPoint(F,V,searchRadius,voxelSize,plotOn);

%%
% Visualize geometry interior/boundary label image
sv3(M,voxelSize);

%%
% Visualize geometry convoluted interior image
sv3(ML,voxelSize); colormap gjet;

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
