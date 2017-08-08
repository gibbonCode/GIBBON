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




 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
