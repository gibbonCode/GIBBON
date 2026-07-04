%% hemiSphereRegionMesh
% Below is a demonstration of the features of the |ind2patch| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,regionIndSub]=hemiSphereRegionMesh(hemiSphereStruct);|

%% Description
% The |hemiSphereRegionMesh| function creates the faces (F), vertices (or
% nodes, V) and the region indices (regeionIndSub) for a hemi-sphere
% according to the input structure hemiSphereStruct. The latter defines the
% sphere radius, the number of refinement steps for the regions and the
% number of refinement steps for the mesh. For more information on the
% refinement see the |geoSphere| and |subTri| functions and associated demo
% files. 
% A complete sphere is first represented as an icosahedron which is then
% refined (subtriangulated) hemiSphereStruct.nRefineRegions times (whereby
% for each iteration each triangle is subdevided into 4 triangles). This
% initial subdevision defines the element regions. The next refinement step
% defines the number of triangles for each region. The field
% hemiSphereStruct.nRefineMesh defines how many times each mesh region is
% iteratively subtriangulated. 

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
fontSize=25;
faceAlpha=1;
lineWidth=1;
markerSize=5;

%% Example:  Creating a hemisphere mesh using the |hemiSphereRegionMesh| function
% Defining hemi-sphere parameters
hemiSphereStruct.sphereRadius=1; %Sphere radius
hemiSphereStruct.nRefineRegions=1; %Number of refinement steps for regions
hemiSphereStruct.nRefineMesh=2; %Number of refinement steps for mesh

% Get hemi-sphere mesh
[F,V,regionInd]=hemiSphereRegionMesh(hemiSphereStruct);

%% 
% Plotting results

%Creating a random color for the each mesh region
cmap=hsv(max(regionInd(:)));
cmap=cmap(randperm(size(cmap,1)),:); %scramble colors

hf=cFigure; hold on; 
gtitle('Half dome showing regions with subtriangulated mesh',fontSize);
gpatch(F,V,regionInd);
colormap(cmap);
axisGeom(gca,fontSize);
camlight headlight; 
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
