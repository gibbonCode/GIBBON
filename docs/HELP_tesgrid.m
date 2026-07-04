%% tesgrid
% Below is a demonstration of the features of the |tesgrid| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=tesgrid(xrange,yrange);|
% |[E,V]=tesgrid(xrange,yrange,zrange);|

%% Description
% This function generates a gridded mesh and outputs a tesselation for that
% mesh. The input is similar to that of the |ndgrid| function in that a
% desired range for the x, y, and z coordinates is given. The output
% consists of the mesh simplices. For 2D input only the x and y ranges are
% specified and the output consists of quadrilateral faces and vertices of
% a 2D mesh. For 3D input the z coordinate range is also specified and the
% output consists of hexahedral elements and vertices. 

%% Examples

%%
% Plot settings

fontSize=20;
faceAlpha1=0.8;

%% Creating a quadrilateral gridded mesh using |tesgrid|
% Specifying coordinate ranges
xrange=linspace(-5,5,20);
yrange=linspace(0,11,11);

%% 
% Using |tesgrid| to create a gridded quadrilateral mesh
[F,V]=tesgrid(xrange,yrange);

%%
% Plotting model
cFigure; hold on;
title('A gridded quadrilateral mesh','FontSize',fontSize);
gpatch(F,V,'gw','k',faceAlpha1); %Visualize the mesh
patchNormPlot(F,V); %Show normal directions
axisGeom(gca,fontSize); 
camlight headlight;
drawnow; 

%% Creating a hexahedral gridded mesh using |tesgrid|
% Specifying coordinate ranges
xrange=linspace(-3,3,6);
yrange=linspace(0,11,4);
zrange=linspace(-pi,pi,6);

%% 
% Using |tesgrid| to create a gridded quadrilateral mesh
[E,V]=tesgrid(xrange,yrange,zrange);
[F]=element2patch(E);

%%
% Plotting model
cFigure; hold on;
title('A gridded hexahedral mesh','FontSize',fontSize);
gpatch(F,V,'gw','k',faceAlpha1); %Visualize the mesh
patchNormPlot(F,V); %Show normal directions
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
