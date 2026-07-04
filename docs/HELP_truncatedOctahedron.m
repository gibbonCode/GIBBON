%% truncatedOctahedron
% Below is a demonstration of the features of the |truncatedOctahedron| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C]=truncatedOctahedron(r)|

%% Description
% This function creates the faces (F) and vertices (V) for a
% truncated octahedron. The output C contains a labelling for the faces
% where 0 denotes hexagons and 1 denotes squares. 

%% 
% Plot settings
fontSize=25;

%% Creating a patch model of a rhombic dodecahedron

r=sqrt(5)/4; %Radii, results in a width of 1

[F,V,C]=truncatedOctahedron(r);

%%
% Plotting results5
cFigure; hold on;
title('A truncated octahedron','FontSize',fontSize);
gpatch(F,V,C);
axisGeom(gca,fontSize);
camlight('headlight'); 
icolorbar;
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
