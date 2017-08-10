%% hexMeshSphere
% Below is a demonstration of the features of the |hexMeshSphere| function

%%

clear; close all; clc;

%%
% Plot settings
figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha1=1;
edgeColor=0.25*ones(1,3);
edgeWidth=2;

%% Creating a hollow hexahedral mesh sphere 
% Creating a solid hexahedral mesh sphere

%Control settings
cPar.sphereRadius=10;
cPar.coreRadius=5;
cPar.numElementsMantel=5; 
cPar.numElementsCore=8; 
cPar.makeHollow=0;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

%%
% Plotting sphere model

hf=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('The hexahedral mesh sphere','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

subplot(1,2,2);
title('Cut-view of the mesh','FontSize',fontSize);

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,~]=element2patch(E(L,:),[],'hex8');
patch('Faces',Fs,'Vertices',V,'FaceColor','b','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

drawnow; 

%% Creating a solid hexahedral mesh sphere 
% Creating a solid hexahedral mesh sphere

%Control settings
cPar.sphereRadius=10;
cPar.coreRadius=5;
cPar.numElementsMantel=5; 
cPar.numElementsCore=8; 
cPar.makeHollow=1;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

%%
% Plotting sphere model

hf=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('The hexahedral mesh sphere','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

subplot(1,2,2);
title('Cut-view of the mesh','FontSize',fontSize);

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,~]=element2patch(E(L,:),[],'hex8');
patch('Faces',Fs,'Vertices',V,'FaceColor','b','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;

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
