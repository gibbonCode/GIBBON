%% hexMeshBox
% Below is a demonstration of the features of the |hexMeshBox| function

%%
clear all; close all; clc; 

%%
% PLOT SETTINGS
fig_color='w'; fig_colordef='white';
fontSize=20;
cmap=gray(250);
falpha=1;
patch_types={'sx','sy','sz','v'};
ptype=3;
no_slices=4;
mark_siz1=25;
mark_siz2=25;
mark_siz3=15;
edgeWidth=2;
edgeColor=0.7*ones(1,3);
faceAlpha1=1;
faceAlpha2=0.5;

%% CREATING A MESHED BOX
boxDim=[5 6 7];
boxEl=[4 5 6];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E=meshStruct.E;
V=meshStruct.V;
F=meshStruct.F;
Fb=meshStruct.Fb;
faceBoundaryMarker=meshStruct.faceBoundaryMarker;

%%
% Plotting model
hf1=figuremax(fig_color,fig_colordef);

subplot(1,2,1);
title('Box boundaries and normals','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fb,V,1);
colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

subplot(1,2,2);
title('Box model','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hpm=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','g','FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
axis equal; view(3); axis tight;  grid on; 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)