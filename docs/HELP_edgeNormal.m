%% edgeNormal
% Below is a demonstration of the features of the |edgeNormal| function

%%
clear; close all; clc;

%% Syntax
% |[NE,VE]=edgeNormal(F,V);|

%% Description
% This function computes the edge normal vectors (NE), and optionally the
% position vectors for these normal vectors (VE), for the input patch data
% defined by faces (F) and vertices (V). 

%% Examples

%% Calculate edge normals for patch data for a single face

%%
% Create example patch surface data

%A single triangle
V=[0 0 0; 1 0 0; 1 1 0];
F=[1 2 3];

%%
% Get edge normals

[NE,VE]=edgeNormal(F,V); %Get edge normals

%%
% Visualize edge normals

%Use mean edge length devided by 4 as the length for visualizing edge
%vectors.
d=mean(patchEdgeLengths(F,V))/4; 

%Create example color data for faces and edges
C=(1:size(F,1))'; %Color data for each face equal to face numbers
CE=[C(:,ones(size(F,2),1))]'; %Get corresponding colors for edges
CE=CE(:);

%Visualize
cFigure; 
gpatch(F,V,C,'k');
hp(1)=patchNormPlot(F,V); %Plot face normals
hp(2)=quiverVec(VE,NE,d,CE,'r'); %Plot edge normals
legend(hp,{'Face normals','Edge normals'});
axisGeom;
camlight headlight; 
colormap viridis
drawnow; 

%% Calculate edge normals for general patch data

%%
% Create example patch surface data

[F,V]=geoSphere(1,1); %A geodesic sphere triangulation

%%
% Get edge normals

[NE,VE]=edgeNormal(F,V); %Get edge normals

%%
% Visualize edge normals

%Use mean edge length devided by 4 as the length for visualizing edge
%vectors.
d=mean(patchEdgeLengths(F,V))/4; 

%Create example color data for faces and edges
C=(1:size(F,1))'; %Color data for each face equal to face numbers
CE=[C(:,ones(size(F,2),1))]'; %Get corresponding colors for edges
CE=CE(:);

%Visualize
cFigure; 
gpatch(F,V,C,'k');
hp(1)=patchNormPlot(F,V); %Plot face normals
hp(2)=quiverVec(VE,NE,d,CE,'r'); %Plot edge normals
legend(hp,{'Face normals','Edge normals'});
axisGeom;
camlight headlight; 
colormap viridis
drawnow; 

%%
