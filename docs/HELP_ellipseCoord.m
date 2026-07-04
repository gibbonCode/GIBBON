%% ellipseCoord
% Below is a demonstration of the features of the |ellipseCoord| function

%%
clear; close all; clc;

%% Syntax
% |[V]=ellipseCoord(A,t);|

%% Description 
% Calculates ellipse coordinates for the angles in t based on the vector A
% which defines the centre coordinates, the radii and the angle
% respectively. 

%% Examples 
% 

%Define example input data

%Angles
t=linspace(0,2*pi);

%Centre coordinates
xc=2;
yc=2; 

%Radii
r1=1;
r2=2;

%Orientation angle
a=0.25*pi; 

%Compose A
A=[xc yc r1 r2 a];

%Compute ellipse coordinates
V=ellipseCoord(A,t);

%%
% Visualize ellipse

cFigure; 
plotV(V,'r.-','MarkerSize',25,'LineWidth',3);
axisGeom; view(2);
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
