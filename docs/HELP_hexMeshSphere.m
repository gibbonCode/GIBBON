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
