%% honeyCombMesh
% Below is a demonstration of the features of the |honeyCombMesh| function

%%
clear; close all; clc; 

%%
% PLOT SETTINGS
figColor='w'; 
figColorDef='white';
fontSize=10;
lineWidth1=3;
faceAlpha1=0.9;

%% 
% Control parameters

%Desired mesh point spacing
pointSpacing=2;

%Mesh region extrema
maxV=[10 10];
minV=[-10 -10];

%% CREATING A HONEY-COMB MESH

[Fh,Vh]=honeyCombMesh(minV,maxV,pointSpacing);

%%
% Plottting model
hf1=figuremax(figColor,figColorDef);
title('The honey-comb mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
C=rand(size(Fh,1),1);
hpm=patch('Faces',Fh,'Vertices',Vh,'EdgeColor','k','FaceColor','flat','CData',C,'FaceAlpha',faceAlpha1,'lineWidth',lineWidth1);
% [hp]=patchNormPlot(Fh,Vh,1);
colormap autumn; 
axis equal; view(3); axis tight;  grid on; 
set(gca,'FontSize',fontSize);
drawnow;
view(2);

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>