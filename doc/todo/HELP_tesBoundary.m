%% tesBoundary
% Below is a demonstration of the features of the |tesBoundary| function

%%
clear; close all; clc;

%% Syntax
% |[indBounary]=tesBoundary(F,V);|

%% Description 
% UNDOCUMENTED 

%% Examples 
% 
%%
% Plot settings
faceAlpha=0.5; 
fontSize=20; 

%% Example 1, Obtaining boundary faces of a hexahedral mesh

boxDim=[5 5 1];
boxEl=[5 5 1];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E=meshStruct.E;
V=meshStruct.V;
C=(1:1:size(E,1))';

E=E([1:3 7],:);

% [V,~]=platonic_solid(2,1); %Create a cube
% E=[1:8]; %
% [E,V]=subHex(E,V,2);
% C=(1:1:size(E,1))';
% 
% % [~,indUni,indFix]=unique(pround(V,5),'rows');
% % V=V(indUni,:);
% % E=indFix(E); 
% 
[F,CF]=element2patch(E,C,'hex8');

%%
% Plotting model
cFigure; hold on;
title('Hexahedral mesh','FontSize',fontSize);

gpatch(F,V,0.5*ones(1,3),'k',faceAlpha);

axisGeom(gca,fontSize); 
camlight headlight;

drawnow;

%%

[indB]=tesBoundary(F,V);

Ft=[F(indB,[1 2 3]); F(indB,[3 4 1]);];

%%
% Plotting model
cFigure; hold on;
title('Hexahedral mesh','FontSize',fontSize);

% gpatch(F(indB,:),V,0.5*ones(1,3),'k',faceAlpha);
gpatch(Ft,V,0.5*ones(1,3),'k',faceAlpha);
% gpatch(F,V,0.5*ones(1,3),'k',faceAlpha);
% plotV(V(indB,:),'r.');

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
