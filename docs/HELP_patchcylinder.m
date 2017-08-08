%% patchcylinder
% Below is a demonstration of the features of the |patchcylinder| function

%% Syntax
% |[F,V]=patchcylinder(r,nr,h,nz,ptype);|

%% Description
% Use |patchcylinder| to generate triangulated spheres with nearly geodesic
% triangle distributions. The density of the triangulation can be
% controlled through a particular choice of n (number of mesh refinement
% steps).

%% Examples

clear; close all; clc;

%% 
% Plot Settings
fontSize=15;
faceAlpha=0.8;
edgeColor=0.*ones(1,3);
edgeWidth=2;

%% Building a patched cylinder model
% Defining patched cylinders with different mesh types

cylRaduis=1; %Cylinder radius
numRadial=12; %Number of elements in the circumferential direction
cylHeight=3; %height
numSteps=7; %Number of elements in the height direction
meshTypes={'quad','tri_slash','tri'}; %Patch Types

pColors=gjet(numel(meshTypes));

%%
% Creating and visualizing patch data

cFigure; 
for q=1:1:numel(meshTypes)
    [F,V]=patchcylinder(cylRaduis,numRadial,cylHeight,numSteps,meshTypes{q}); 
    subplot(1,3,q); hold on;
    title([meshTypes{q},' type cylinder'],'FontSize',fontSize,'Interpreter','none');
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    gpatch(F,V,pColors(q,:));
    patchNormPlot(F,V);    
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  grid on;
    camlight headlight;
end
drawnow;

%% Using an input structure instead

%%
% Creating input structure
inputStruct.cylRadius=1;
inputStruct.numRadial=15;
inputStruct.cylHeight=3;
inputStruct.numHeight=11;
inputStruct.meshType='tri';

%%
% Derive patch data for a cylinder
[F,V]=patchcylinder(inputStruct); 

%%
% Visualizing cylinder model

cFigure; 
hold on;
gpatch(F,V,'g');
% patchNormPlot(F,V);  
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Creating a closed cylinder model

%%
% Indices for the top and bottom points can be obtained as follows
indTop=inputStruct.numHeight:inputStruct.numHeight:size(V,1);
indBottom=1:inputStruct.numHeight:size(V,1);

%%
% The top and bottom can be meshed using |regionTriMesh2D|

[Ft,Vt]=regionTriMesh2D({V(indTop,[1 2])},[],0);
Vt(:,3)=mean(V(indTop,3));

[Fb,Vb]=regionTriMesh2D({V(indBottom,[1 2])},[],0);
Vb(:,3)=mean(V(indBottom,3));
Fb=fliplr(Fb);

%%
% Visualizing meshed regions

%Curves
plotV(V(indTop,:),'b.-','markerSize',25,'lineWidth',2);
plotV(V(indBottom,:),'r.-','markerSize',25,'lineWidth',2);

%Caps
gpatch(Ft,Vt,'r');
% patchNormPlot(Ft,Vt);  
gpatch(Fb,Vb,'b');
% patchNormPlot(Fb,Vb);  
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
