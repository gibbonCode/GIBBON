%% axisGeom
% Below is a demonstration of the features of the |axisGeom| function

%%
clear; close all; clc;

%% Syntax
% |axisGeom(h,fontSize);|

%% Description 
% This function sets axis properties commonly used for viewing geometry in
% 3D:  
% 
% |axis equal; axis vis3d; axis tight;|   
% 
% The optional inputs include the axis handle and the desired font size. 

%% Examples 
% 

%% Setting axis parameters for 3D geometry viewing

%%
% Create something to visualize
[F,V]=graphicsModels(9);

%%
% Visualize 

cFigure; hold on; 
gpatch(F,V,'w','k');

axisGeom; %Set axis properties

camlight headlight;
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
