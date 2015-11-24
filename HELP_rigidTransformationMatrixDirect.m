%% HELP_rigidTransformationMatrixDirect
% Below is a demonstration of the features of the |rigidTransformationMatrixDirect| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceAlpha=1;
edgeColor=0.3*ones(1,3);
edgeWidth=1.5;

%% 
% Load example patch data
[F,V1]=parasaurolophus;

V2=V1;
% V2(:,1)=V2(:,1)/2; 

%%
% Create a test tranformation matrix

%Define a rotation
[R,~]=euler2DCM([0.25*pi 0.75*pi 0.1*pi]); 

%Build a tranformation matrix
T=R; 

%Add translation
T(:,4)=[-0.1 2.1 0.5]; 
T(4,:)=0; 
T(4,4)=1;

%Transform
VV=V2; 
VV(:,4)=1; 
VT=(T*VV')'; 
V2=VT(:,[1 2 3]);

%%
% Plotting data

hf=figuremax(fig_color,fig_colordef);
title('The original (green) and transformed surface (red)','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F,'Vertices',V1,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);

hp=patch('Faces',F,'Vertices',V2,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);

camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 
drawnow; 

%% Get the transformation matrix for the point matched data using |rigidTransformationMatrixDirect|

[Tf]=rigidTransformationMatrixDirect(V1,V2);

%Compare T and Tf
disp(T);
disp(Tf);

%The residuals
r=T-Tf;
disp(r);

%Inverse map coordinates
V1f=(Tf\VT')'; 
V1f=V1f(:,[1 2 3]);

%%
% Plotting data

hf=figuremax(fig_color,fig_colordef);
title('The original (green) and inverse mapped surface (red)','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F,'Vertices',V1,'FaceColor','g','FaceAlpha',0.5,'lineWidth',edgeWidth,'edgeColor',edgeColor);

hp=patch('Faces',F,'Vertices',V1f,'FaceColor','none','FaceAlpha',0.5,'lineWidth',edgeWidth*2,'edgeColor','r');

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 
drawnow; 

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>