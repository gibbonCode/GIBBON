%% quadPlate
% Below is a demonstration of the features of the |quadPlate| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=quadPlate(plateDim,plateEl)|

%% Description
% This function generates a quadrilateral mesh for a rectangular place

%% Examples

%%
% Plot settings
fontSize=20;

%% Example: Creating a meshed plate

%Define input parameters
plateDim = [2 2]; %Plate dimensions in X and Y direction
plateEl  = [6 7]; %Number of elements in X and Y direction

%Use |quadPlate| to create a quadrilateral mesh
[F,V]=quadPlate(plateDim,plateEl);

%%
% Plotting model
cFigure; hold on;
gpatch(F,V,'bw','k');
% patchNormPlot(F,V);
axisGeom(gca,fontSize); 
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
