%% HELP_euler2DCM
% Below is a demonstration of the features of the |euler2DCM| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceAlpha=1;
edgeColor=0*ones(1,3);
edgeWidth=1;

%% Creating a rotation matrix using Euler angles
% Search wikipedia for information on Euler angles to define a rotation or
% direction cosine matrix (DCM). 

%% 
% Load example patch data
[F,V]=parasaurolophus;

%%
% Defining sets of Euler angles for X, Y and Z axis rotation

E1=[0.5*pi 0 0]; 
E2=[0 0.5*pi 0]; 
E3=[0 0 0.5*pi]; 
E4=[0.25*pi 0.25*pi 0.25*pi]; 

%%
% Use |euler2DCM| function to define the rotation matrices

[R1]=euler2DCM(E1);
[R2]=euler2DCM(E2);
[R3]=euler2DCM(E3);
[R4]=euler2DCM(E4);

%%
% Rotate the coordinates

V1=(R1*V')'; 
V2=(R2*V')'; 
V3=(R3*V')'; 
V4=(R4*V')'; 

%%
% Plotting data

hf=figuremax(fig_color,fig_colordef);

subplot(2,2,1);
title('X-axis rotation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',F,'Vertices',V1,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 

subplot(2,2,2);
title('Y-axis rotation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',F,'Vertices',V2,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 

subplot(2,2,3);
title('Z-axis rotation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',F,'Vertices',V3,'FaceColor','b','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 

subplot(2,2,4);
title('Off-axis rotation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',F,'Vertices',V4,'FaceColor','y','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 

drawnow; 

%%
% A second output can also be requested which is the inverse rotation matrix. 
[R,Ri]=euler2DCM([randn(1,3)*pi]);
R
Ri

%%
% i.e. such that the following :

Vr=(R*V')'; %The rotated coordinates
Vn=(Ri*Vr')'; %The normal coordinates after transforming back the rotated coordinates using inverse matrix

%%
% Note that the sum of squared differences for instance is nearly zero 
D=sum((V(:)-Vn(:)).^2)

%% Creating multiple rotation matrices
% It is possible to define multiple rotation matrices at once by specifying
% a multi-row angle set

E=[0.25*pi 0 0; 0 0.5*pi 0]; %E.g. two angle sets are specified, 1 for each row

%%
% In this case the rotation matrices are stacked in the 3rd dimension

[R]=euler2DCM(E)

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
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
