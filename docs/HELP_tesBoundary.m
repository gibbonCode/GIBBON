%% tesBoundary
% Below is a demonstration of the features of the |tesBoundary| function

%%
clear; close all; clc;

%% Syntax
% |[indBoundary]=tesBoundary(F);|

%% Description 
% This function obtains the indices of the boundary faces based on the
% input faces F. Boundary faces are those that occur only once (in terms of
% nodes involved, e.g. [1 2 3 4] is deemed the same face as [2 3 4 1]. 

%% Examples 
% 

%%
% Plot settings

fontSize=20;
faceAlpha1=0.5;

%% Example 1: Get hex mesh boundary faces

%% 
% Example geometry

boxDim=[6 4 4];
boxEl=[5 3 3];

[meshStruct]=hexMeshBox(boxDim,boxEl,2);

E=meshStruct.elements; %Elements
V=meshStruct.nodes; %Nodes
F=meshStruct.faces; %Mesh faces
%The hexMeshBox function provides boundary faces already but tesboundary is used in this example
% Fb=meshStruct.facesBoundary; %Boundary faces

%% 
% Get boundary faces

indBoundary=tesBoundary(F);

Fb=F(indBoundary,:);

%%
% Plotting model
cFigure;
subplot(1,2,1); hold on;
title('All faces','FontSize',fontSize);
gpatch(F,V,'w','k',faceAlpha1,3);
axisGeom(gca,fontSize); camlight headlight; 

subplot(1,2,2); hold on;
title('Boundaries faces','FontSize',fontSize);
gpatch(Fb,V,'w','k',faceAlpha1,3);
axisGeom(gca,fontSize); camlight headlight; 

drawnow; 

%% Example 1: Get tet mesh boundary faces

%% 
% Example geometry

boxDim=[6 4 4];
pointSpacing=2;

[meshStruct]=tetMeshBox(boxDim,pointSpacing);

E=meshStruct.elements; %Elements
V=meshStruct.nodes; %Nodes
F=meshStruct.faces; %Mesh faces
%The tetMeshBox function provides boundary faces already but tesboundary is used in this example
% Fb=meshStruct.facesBoundary; %Boundary faces

%% 
% Get boundary faces

indBoundary=tesBoundary(F);

Fb=F(indBoundary,:);

%%
% Plotting model
cFigure;
subplot(1,2,1); hold on;
title('All faces','FontSize',fontSize);
gpatch(F,V,'w','k',faceAlpha1,3);
axisGeom(gca,fontSize); camlight headlight; 

subplot(1,2,2); hold on;
title('Boundaries faces','FontSize',fontSize);
gpatch(Fb,V,'w','k',faceAlpha1,3);
axisGeom(gca,fontSize); camlight headlight; 

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
