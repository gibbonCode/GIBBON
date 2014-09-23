%% honeyCombMesh
% Below is a demonstration of the features of the |honeyCombMesh| function

%%
clear all; close all; clc; 

%%
% PLOT SETTINGS
fig_color='w'; fig_colordef='white';
font_size=20;
cmap=gray(250);
falpha=1;
patch_types={'sx','sy','sz','v'};
ptype=3;
no_slices=4;
mark_siz1=25;
mark_siz2=25;
mark_siz3=15;
line_width1=3;
F_alpha1=0.9;
F_alpha2=0.8;

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
hf1=figuremax(fig_color,fig_colordef);
title('The honey-comb mesh','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size); zlabel('Z','FontSize',font_size);
hold on;
C=rand(size(Fh,1),1);
hpm=patch('Faces',Fh,'Vertices',Vh,'EdgeColor','k','FaceColor','flat','CData',C,'FaceAlpha',F_alpha1,'lineWidth',line_width1);
% [hp]=patchNormPlot(Fh,Vh,1);
colormap autumn; 
axis equal; view(3); axis tight;  grid on; 
set(gca,'FontSize',font_size);
drawnow;
view(2);
%% 
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)