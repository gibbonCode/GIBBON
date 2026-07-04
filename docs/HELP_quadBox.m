%% quadBox
% Below is a demonstration of the features of the |quadBox| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,faceBoundaryMarker]=quadBox(boxDim,boxEl);|

%% Description
% This function generates the patch data for a quadrilateral mesh
% rectangular box. The patch data consists of the faces (F), the vertices
% (V), and face color/label data (faceBoundaryMarker). The input parameters
% boxDim (1x3 vector) and BoxEl (1x3 vector) define the dimensions of the
% box and the number of elements to use in all 3 directions. 

%% Examples

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;

%% Creating a quadrilateral mesh of a box

%% 
% Specifying dimensions and number of elements for each direction
boxDim=[4 5 6]; %Size of the box in each direction
boxEl=[3 4 5]; %Number of elements per direction 

%%
% Using |quadBox| to build the patch model

[F,V,faceBoundaryMarker]=quadBox(boxDim,boxEl);

%%
% Plotting model
cFigure; hold on;
title('Box quadrilateral faces and normals','FontSize',fontSize);

gpatch(F,V,faceBoundaryMarker,'k',faceAlpha1);
% patchNormPlot(F,V);

colormap(gjet(6)); icolorbar; 
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
