%% tetMeshBox
% Below is a demonstration of the features of the |tetMeshBox| function

%%
clear; close all; clc;

%% Syntax
% |[meshStruct]=tetMeshBox(boxDim,pointSpacing);|

%% Description
% This function generates a mesh structure containing element and node data
% for a tetrahedral element meshed box. The box dimensions in each of the 3
% directions are based on the boxDim input (1x3 vector). The number of
% elements in each direction is based on the pointSpacing input.

%% Examples

%%
% Plot settings

fontSize=20;
faceAlpha1=0.8;

%% CREATING A MESHED BOX
boxDim=[5 6 7]; % Box dimenstions
pointSpacing=1; 

[meshStruct]=tetMeshBox(boxDim,pointSpacing);

%%
% Acces output fields
E=meshStruct.elements;
V=meshStruct.nodes;
F=meshStruct.faces;
Fb=meshStruct.facesBoundary;
faceBoundaryMarker=meshStruct.boundaryMarker;

%%
% Plotting model
cFigure;
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha1);
% patchNormPlot(Fb,V);

axisGeom(gca,fontSize); 
colormap(gjet(6)); icolorbar; 
drawnow; 

%%
% Visualizing model internal mesh with |meshView|

meshView(meshStruct); 

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
