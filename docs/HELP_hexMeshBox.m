%% hexMeshBox
% Below is a demonstration of the features of the |hexMeshBox| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=20;
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
cFigure;
subplot(1,2,1);
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha1);
% [hp]=patchNormPlot(Fb,V,1); %Display face normals

axisGeom(gca,fontSize); 
colormap(gjet(6)); icolorbar; 


subplot(1,2,2);
title('Box hexahedral mesh','FontSize',fontSize);
hold on;

gpatch(F,V,0.5*ones(1,3),'k',faceAlpha2);
patchNormPlot(F,V);
axisGeom(gca,fontSize); 
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