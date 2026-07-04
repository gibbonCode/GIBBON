%% triIncenter
% Below is a demonstration of the features of the |triIncenter| function

%%
clear; close all; clc;

%% Syntax
% |[P]=triIncenter(F,V);|

%% Description 
% This function computes the incentres P for the input triangulation
% definced by the faces F and vertices V. 

%% Examples 
% 

%%
% Plot settings for examples
markerSize=15; 

%% Example 1: Computing incenters for an exampe triangulated mesh
%

% Get vertices and faces for example geometry
[F,V]=stanford_bunny;

%%
% Use |triIncenter| to compute the incentres

[P]=triIncenter(F,V);

%%

cFigure; hold on; 
gpatch(F,V,'w','o',1,1);
hp=plotV(P,'r.','MarkerSize',markerSize);
legend(hp,'Triangle incenters');
axisGeom; camlight headlight;
gdrawnow; 

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
