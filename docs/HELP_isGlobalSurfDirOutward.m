%% isGlobalSurfDirOutward
% Below is a demonstration of the features of the |isGlobalSurfDirOutward| function

%%
clear; close all; clc;

%% Syntax
% |[L]=isGlobalSurfDirOutward(F,V);|

%% Description 
% This function returns a boolean denoting wether surface normals are
% pointing outward (1) (which would result in a positive volume being computed)
% or inward (2) (which would result in a negative volume being computed). 
% The method assumes that the normal directions are coherent across the
% surface. 

%% Examples 
% 

%%
% Example surface
[F,V]=stanford_bunny;

%% 
% Computing a logic denoting whether the normals point outward or inwards.
% For the example surface the normals point outward, so the logical should
% return a true

% This one should be true
[L]=isGlobalSurfDirOutward(F,V)

%%
% However inverting the surface causes the normal directions to point
% inward, hence a false is returned. 

% Reversed should be false
[L]=isGlobalSurfDirOutward(fliplr(F),V)

%%

cFigure; 
gpatch(F,V);
patchNormPlot(F,V);
axisGeom; camlight headlight; 
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
