%% rhombicDodecahedron
% Below is a demonstration of the features of the |rhombicDodecahedron| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=rhombicDodecahedron(r);|

%% Description 
% This function creates the faces (F) and vertices (V) for a
% rhombic-dodecahedron. 

%% Examples 
% 

%% 
% Plot settings
fontSize=15;
faceAlpha1=0.3;
edgeColor='k';
lineWidth1=3;
markerSize=55;
faceColor=0.5*ones(1,3);

%% Creating a patch model of a rhombic dodecahedron

r=sqrt(2)/2; %Radii, results in a width of 1

[F,V]=rhombicDodecahedron(r);

R=euler2DCM([0 0 0.25*pi]);
V=V*R;

%%
% Plotting results
cFigure;
title('A rhombic dodecahedron','FontSize',fontSize);
hold on;

gpatch(F,V,faceColor,edgeColor,faceAlpha1,lineWidth1);
plotV(V,'k.','MarkerSize',markerSize);
patchAnnotate(F,V,[],'fontSize',fontSize);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow;

%%
% Plotting results
cFigure; hold on;
gpatch(F,V,'bw',edgeColor,1,lineWidth1);
axisGeom(gca,fontSize);
camlight('headlight'); 
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
