%% hemiSphereCylMesh
% Below is a demonstration of the features of the |hemiSphereCylMesh| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=hemiSphereCylMesh(inputStruct);|

%% Description
% This function generates patch data for a cylinder which ends with a
% hemispherical head. 

%% Examples

%%
% Plot settings
figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha=1;
lineWidth=1;
markerSize1=10;

%% 
% Creating 3 example input structures 

inputStruct1.sphereRadius=3; %Sphere radius
inputStruct1.nRefine=2; %Number of refinement steps for sphere regions
inputStruct1.cylinderHeight=5; %Height of the cylindrical part
inputStruct1.cylinderStepSize=[]; %Aproximate desired node spacing for cylindrical part, empty uses spacing of hemi-sphere mesh
inputStruct1.patchType='tri';

inputStruct2.sphereRadius=3; %Sphere radius
inputStruct2.nRefine=2; %Number of refinement steps for sphere regions
inputStruct2.cylinderHeight=5; %Height of the cylindrical part
inputStruct2.cylinderStepSize=1; %Aproximate desired node spacing for cylindrical part, empty uses spacing of hemi-sphere mesh
inputStruct2.patchType='tri_slash';

inputStruct3.sphereRadius=3; %Sphere radius
inputStruct3.nRefine=2; %Number of refinement steps for sphere regions
inputStruct3.cylinderHeight=5; %Height of the cylindrical part
inputStruct3.cylinderStepSize=[]; %Aproximate desired node spacing for cylindrical part, empty uses spacing of hemi-sphere mesh
inputStruct3.patchType='quad';

%% CREATING A SURFACE TRIANGULATION COMPOSED OF A MERGED HEMI-SPHERE AND CYLINDER
[F1,V1]=hemiSphereCylMesh(inputStruct1);
[F2,V2]=hemiSphereCylMesh(inputStruct2);
[F3,V3]=hemiSphereCylMesh(inputStruct3);

%% PLOTTING MODEL

cFigure; 
subplot(1,3,1); hold on; 
title('auto point spacing, tri','FontSize',fontSize);
gpatch(F1,V1,'gw');
patchNormPlot(F1,V1);
axisGeom(gca,fontSize); 
camlight headlight; 

subplot(1,3,2); hold on; 
title('custom point spacing, tri_slash','FontSize',fontSize,'interpreter','none');
gpatch(F2,V2,'gw');
patchNormPlot(F2,V2);
axisGeom(gca,fontSize); 
camlight headlight; 

subplot(1,3,3); hold on; 
title('auto point spacing, quad','FontSize',fontSize);
gpatch(F3,V3,'gw');
patchNormPlot(F3,V3);
axisGeom(gca,fontSize); 
camlight headlight; 

drawnow; 

%% Example: Face color data
% Colors can be used as handles to surface components. The option color
% data output can be used to seperate the hemisphere form the cylindrical
% part of the mesh. 

%%
% Create mesh
[F1,V1,C1]=hemiSphereCylMesh(inputStruct1);

%%
% Visualize color data

cFigure; hold on; 
title('Face color data','FontSize',fontSize);
gpatch(F1,V1,C1);
axisGeom(gca,fontSize);
colormap gjet; icolorbar; 
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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
