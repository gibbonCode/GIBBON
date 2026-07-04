%% triBox
% Below is a demonstration of the features of the |triBox| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=15;
faceAlpha1=0.5;

%% Creating a triangulated mesh of a box

%% 
% Specifying dimensions and number of elements for each direction
boxDim=[4 5 6]; %Width in each direction
pointSpacing=1; %Desired point spacing

%%
% Using |triBox| to build the patch model

[F,V,faceBoundaryMarker]=triBox(boxDim,pointSpacing);

%%
% Visualisation
cFigure; hold on;
title('Box triangular faces and normals','FontSize',fontSize);

gpatch(F,V,faceBoundaryMarker,'k',1);
patchNormPlot(F,V);

axisGeom(gca,fontSize);
colormap(gjet(6)); icolorbar; 
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
