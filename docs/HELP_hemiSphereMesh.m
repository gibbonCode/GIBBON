%% hemiSphereMesh
% Below is a demonstration of the features of the |hemiSphereMesh| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C]=hemiSphereMesh(nRefineSteps,sphereRadius,closeOpt);|

%% Description 
% Creates the patch data (faces F, vertices V, and colordata C) for a
% hemisphere with a radius equal to sphereRadius. The hemisphere is
% refined nRefineSteps times and is closed if closeOpt==1
%
%% Examples 
% 

nRefineSteps=2; %Number of refinement steps
sphereRadius=1; %Radius
closeOpt=1; %Option to close bottom of hemisphere

[F,V,C]=hemiSphereMesh(nRefineSteps,sphereRadius,closeOpt); %Construct hemi-shere mesh

%%
% Visualize mesh

cFigure; 
subplot(1,2,1); hold on;
gpatch(F,V,'w','k',1,2);
axisGeom;
camlight headlight;

subplot(1,2,2); hold on;
gpatch(F,V,C,'none',0.5);
axisGeom;
camlight headlight;
colormap viridis; icolorbar; 

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
