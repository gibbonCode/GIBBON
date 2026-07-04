%% rhombicDodecahedronHex
% Below is a demonstration of the features of the |rhombicDodecahedronHex| function

%%
clear; close all; clc;

%% Syntax
% |[E,V]=rhombicDodecahedronHex(r);|

%% Description 
% This function generates a rhombic dodecahedron domain which is subdevided
% into hexahedral elements. The input is the bounding sphere radius, the
% outputs are the hexahedral element array and the vertices. 

%% Examples 
% 

%% 
% Plot settings
fontSize=15;
faceAlpha1=0.5;
edgeColor='k';
lineWidth1=3;
markerSize=55;

%% Creating a patch model of a rhombic dodecahedron

r=sqrt(2)/2; %Radii, results in a width of 1

[E,V]=rhombicDodecahedronHex(r);
C=(1:1:size(E,1))';
[F,CF]=element2patch(E,C); 

%%
% Plotting results
cFigure;
title('A hexahedral mesh of a rhombic dodecahedron','FontSize',fontSize);
hold on;

gpatch(F,V,CF,'k',faceAlpha1,lineWidth1);
plotV(V,'k.','MarkerSize',markerSize);

axisGeom(gca,fontSize);
view(-10,25);
camlight('headlight'); 
colormap gjet; icolorbar;
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
