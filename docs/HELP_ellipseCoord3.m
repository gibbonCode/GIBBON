%% ellipseCoord3
% Below is a demonstration of the features of the |ellipseCoord3| function

%%
clear; close all; clc;

%% Syntax
% |[V]=ellipseCoord3(e,t);|

%% Description 
% Calculates ellipse coordinates for the angles in t based on the input
% structure e which contains the following fields: 
% radii, a 2x1 array
% axes, a 3x3 rotation matrix
% centre the ellipse centre coordinates

%% Examples 
% 

%Define example input data

%Angles
t=linspace(0,2*pi);

%Centre coordinates
Vc=[2 2 2];

%Radii
r=[1 2];

%Rotation matrix
Q=euler2DCM([0 -0.25*pi 0.25*pi]); 

%Compose input structure
e.centre=Vc; 
e.radii=r;
e.axes=Q;

%Compute ellipse coordinates
V=ellipseCoord3(e,t);

%%
% Visualize ellipse

cFigure; 
plotV(V,'r.-','MarkerSize',25,'LineWidth',3);
axisGeom;
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
