%% subTet
% Below is a demonstration of the features of the |subTet| function

%% Syntax
% |[Es,Vs]=subTet(E,V,splitMethod);|

%% Description 
% 
%% Examples 

clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceColor1='g';
faceColor2='r';
faceAlpha1=0.3;
faceAlpha2=1;
edgeColor=0.*ones(1,3);
edgeWidth=2;
markerSize=2;
cMap=gjet(250);

%% 
% Create a test tesselation
r=1; %Radius of tetrahedron circumsphere
[V,~]=platonic_solid(1,r);
E=[1:4];

%% Example creating sub-tetrahedrons

%Method 1
[E1,V1]=subTet(E,V,1);
 
%Method 2
[E2,V2]=subTet(E,V,2);

%% Visualization

C1=(1:size(E1,1))'; %Element colors
[F1,CF1]=element2patch(E1,C1); %Patch data for plotting

C2=(1:size(E2,1))'; %Element colors
[F2,CF2]=element2patch(E2,C2);  %Patch data for plotting

cFigure;
subplot(1,2,1); 
title('Method 1','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F1,'Vertices',V1,'FaceColor','flat','CData',CF1,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
set(gca,'FontSize',fontSize);

subplot(1,2,2); 
title('Method 2','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F2,'Vertices',V2,'FaceColor','flat','CData',CF2,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
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
