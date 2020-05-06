%% hexVol
% Below is a demonstration of the features of the |hexVol| function

%%
clear; close all; clc;

%% Syntax
% |[VE,L]=hexVol(E,V);|

%% Description
% This function computes hexahedral element volumes. The input is the
% element description (E) and the nodes (V). The output is the element
% volumes (always positive) and a logic denoting wheter the element appears
% to be valid (1) or inverted (0). 

%% Example: Computing the volume of hexahedral elements

%%
% Create example geometry

%Creating a single hexahedron
X=[-1;  1; 1; -1; -1;  1; 1; -1;];
Y=[-1; -1; 1;  1; -1; -1; 1;  1;];
Z=[-1; -1;-1; -1;  1;  1; 1;  1;];
Vh=[X(:) Y(:) Z(:)];
Eh=1:8; %The hexahedral element

%Subdevided into 8 smaller elements
[E,V]=subHex(Eh,Vh,1);

%Create faces data for plotting
[F,C]=element2patch(E);

%%
% Computing the volume 
[VE,logicValid]=hexVol(E,V)

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
E_inverted(1,:)=E_inverted(1,[5:8 1:4]);
[VE,logicValid]=hexVol(E_inverted,V)

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
