%% discQuadMesh
% Below is a demonstration of the features of the |discQuadMesh| function

%% Syntax
% |[varargout]=discQuadMesh(nElements,r,f);|

%% Description 
% 
%% Examples 
% 
clear; close all; clc; 

%% PLOT SETTINGS

fontSize=20;

%% Building a quadrilateral circular mesh

ne=12; %Elements in radius
r=2; %Outer radius of disc
f=0.5; %Fraction (with respect to outer radius) where central square appears

%Create the mesh
[Fc,Vc]=discQuadMesh(ne,r,f);

%%
%Visualizing mesh

cFigure; hold on;
title('A circle meshed with quadrilateral elements','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

patch('Faces',Fc,'Vertices',Vc,'EdgeColor','k','FaceColor','g','FaceAlpha',1);

axis equal; axis tight; view(2); box on; grid on; 
set(gca,'FontSize',fontSize);
drawnow;

%% Using other outputs
% The additional outputs Cc and indEdge provide a color labelling for the
% central square and other elements, and the indices for the outer edge
% respectively. 

ne=18; %Elements in radius
r=1; %Outer radius of disc
f=0.6; %Fraction (with respect to outer radius) where central square appears

%Create the mesh
[Fc,Vc,Cc,indEdge]=discQuadMesh(ne,r,f);

%%
%Visualizing mesh

cFigure; hold on;
title('A circle meshed with quadrilateral elements','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

patch('Faces',Fc,'Vertices',Vc,'FaceColor','flat','CData',Cc,'EdgeColor','k','FaceAlpha',1);
plotV(Vc(indEdge,:),'g.-','MarkerSize',25,'LineWidth',2);
colormap(gjet(2)); 
axis equal; axis tight; view(2); box on; grid on; 
set(gca,'FontSize',fontSize);
drawnow;

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>




 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
