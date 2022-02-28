%% tetVol
% Below is a demonstration of the features of the |tetVol| function

%%
clear; close all; clc;

%% Syntax
% |[VE,L]=tetVol(E,V);|

%% Description
% This function computes tetrahedral element volumes. The input is the
% element description (E) and the nodes (V). The output is the element
% volumes (always positive) and a logic denoting wheter the element appears
% to be valid (1) or inverted (0). 

%% Examples

%%
% Plot settings
cMap=gjet(250);
faceAlpha1=1;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';
fontSize=15; 

%% Example: Computing the volume of tetrahedral elements

%%
% Create example geometry

%Creating a hexahedron
X=[-1;  1; 1; -1; -1;  1; 1; -1;];
Y=[-1; -1; 1;  1; -1; -1; 1;  1;];
Z=[-1; -1;-1; -1;  1;  1; 1;  1;];
Vh=[X(:) Y(:) Z(:)];
Eh=1:8; %The hexahedral element

%Convert to a tetrahedron.
[E,V]=hex2tet(Eh,Vh,[],2);

%Create faces data for plotting
[F,C]=element2patch(E);

%%
% Computing the volume 
[VE,logicValid]=tetVol(E,V)

%%
% The summed volume should be 8 for this cube
sum(VE)

%% Visualize mesh and face normals

cFigure; hold on; 
gpatch(F,V,'kw','k',0.5);
patchNormPlot(F,V);
axisGeom;
camlight headlight; 
drawnow; 

%%
% Volumes are always positive but inverted elements have a 0 in the
% inverted logic. In the example below the first element is inverted which
% changes the logic to return a 0 for this element. 

E_inverted=E; 
E_inverted(1,:)=E_inverted(1,[4 1 2 3]);
[VE,logicValid]=tetVol(E_inverted,V)

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
