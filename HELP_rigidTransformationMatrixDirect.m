%% HELP_rigidTransformationMatrixDirect
% Below is a demonstration of the features of the |rigidTransformationMatrixDirect| function

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha=1;
edgeColor='k';
%% 
% Load example patch data
[F,Vd]=parasaurolophus;
 
%Define a rotation
a=[-0.25*pi 0.75*pi 0.1*pi]; 
[R,~]=euler2DCM(a);

V1=Vd+0.5;

V2=tform(R,Vd)-1;

%%
% Plotting data

hf=cFigure;
title('The untransformed surfaces','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F,'Vertices',V1,'FaceColor','g','FaceAlpha',faceAlpha);

hp=patch('Faces',F,'Vertices',V2,'FaceColor','r','FaceAlpha',faceAlpha,'edgeColor',edgeColor);

camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal; 
drawnow; 

%% Get the transformation matrix for the point matched data using |rigidTransformationMatrixDirect|

[M]=rigidTransformationMatrixDirect(V1,V2);

V1f=tform((M),V1);

%%
% Plotting data

hf=cFigure;
title('The green surfaces transformed towards the red','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F,'Vertices',V2,'FaceColor','r','FaceAlpha',0.5,'edgeColor',edgeColor);

hp=patch('Faces',F,'Vertices',V1f,'FaceColor','none','FaceAlpha',0.5,'edgeColor','g');
% hp=patch('Faces',F,'Vertices',V2ff,'FaceColor','none','FaceAlpha',0.5*2,'edgeColor','b');

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