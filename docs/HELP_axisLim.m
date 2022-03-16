%% axisLim
% Below is a demonstration of the features of the |axisLim| function

%%
clear; close all; clc;

%% Syntax
% |axLim=axisLim(V);|

%% Description 
% This function computes appropiate axis limits for the input vertices V.
% The vertices may be an k x l x m array, where by by k is the number of
% vertices, l is the number of dimensions (e.g. 2 or 3), and m is for
% instance a time (or other) dimension. The function returs axLim which are
% appropriate axis limits such that the coordinates in V can be displayed
% appropriately (and in a tight fashion). 

%% Examples

%%
% Create example data

[X,Y,Z]=peaks(25);
[F,V]=surf2patch(X,Y,Z); %Get faces and vertices

%%
% Compute appropriate limits
axLim=axisLim(V) %Axis limits for vertices in V

%%
% Assign axis limits directly

cFigure; 
surf(X,Y,Z);
axisGeom; camlight headlight; 
axis(axisLim(V));
gdrawnow; 

%% Examples 
% 
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
