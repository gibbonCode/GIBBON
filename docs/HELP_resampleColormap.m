%% resampleColormap
% Below is a demonstration of the features of the |resampleColormap| function

%%
clear; close all; clc;

%% Syntax
% |[cmap_i]=resampleColormap(cmap,n);|

%% Description
% The function |resampleColormap| function resamples the input colormap
% cmap using n steps. Resampling is based on linear interpolation. 

%% Examples

%%
% Plot settings
fontSize=15;

%%
% Create example data for visualizations
[X,Y,Z]=peaks(250);

%% Example: Create and resample your own colormap
% Define your own colormap using RGB values (you can search online fpr
% color charts with rbg values, if given as unit8 number, normalize by
% deviding by 255). 

%%
% Create discrete custom colormap
white_rgb=[1 1 1];
red_rgb=[1 0 0];
green_rgb=[0 1 0];
blue_rgb=[0 0 1];
black_rgb=[0 0 0];
rgbData=[white_rgb; red_rgb; green_rgb; blue_rgb; black_rgb];

%%
% Resample colormap e.g. using more intermediate colors

n=250; %Number of color levels
[my_colormap]=resampleColormap(rgbData,n);

%%
% Test colormap

cFigure; hold on; 
title('Custom colormap');
surf(X,Y,Z); shading interp; 
colormap(my_colormap); colorbar;
view(2); axis equal; axis tight; grid on; box on; 
set(gca,'FontSize',fontSize);
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
