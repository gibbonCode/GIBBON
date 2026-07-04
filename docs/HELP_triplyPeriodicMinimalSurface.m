%% triplyPeriodicMinimalSurface
% Below is a demonstration of the features of the |triplyPeriodicMinimalSurface| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C]=triplyPeriodicMinimalSurface(inputStruct);|

%% Description
% This function generates Stochastic Bicontinuous Microstructures
%
% Input structure and default values:
%
% inputStruct.L=1; % characteristic length
% inputStruct.Ns=80; % number of sampling points
% inputStruct.anisotropyFactors=[1 1 1]; 
% inputStruct.isocap=1; %Option to cap the isosurface
% inputStruct.surfaceCase='g'; %Surface type
% inputStruct.numPeriods=[1 1 1]; %Number of periods in each direction
% inputStruct.levelset=0.5; %Isosurface level


%% Examples

%%
% Plot settings
cMap=parula(250);
faceAlpha1=1;
faceAlpha2=0.5;
edgeColor1='none';
edgeColor2='none';
fontSize=15; 

%% Example 1: A closed surface
inputStruct.L=1; % characteristic length
inputStruct.Ns=80; % number of sampling points
inputStruct.anisotropyFactors=[1 1 1]; 
inputStruct.isocap=1; %Option to cap the isosurface
inputStruct.surfaceCase='g'; %Surface type
inputStruct.numPeriods=[3 3 3]; %Number of periods in each direction
inputStruct.levelset=0; %Isosurface level
inputStruct.surfaceSide=-1; % 0=both, 1="above" 1, -1="below"

%% 
% Create triply periodic minimal surface
[F,V,C]=triplyPeriodicMinimalSurface(inputStruct);

%%
% Using grouping to keep only largest group
groupOptStruct.outputType='label';
[G,~,groupSize]=tesgroup(F,groupOptStruct); %Group connected faces
[~,indKeep]=max(groupSize); %Index of largest group

%Keep only largest group
F=F(G==indKeep,:); %Trim faces
C=C(G==indKeep,:); %Trim color data 
[F,V]=patchCleanUnused(F,V); %Remove unused nodes

%%
% Visualize surface

cFigure; 
gpatch(F,V,C,'none');
% patchNormPlot(F,V);
axisGeom; camlight headlight; 
colormap gjet; icolorbar;
gdrawnow;

%% Example 2: A non-closed surface
inputStruct.L=1; % characteristic length
inputStruct.Ns=80; % number of sampling points
inputStruct.anisotropyFactors=[1 1 1]; 
inputStruct.isocap=0; %Option to cap the isosurface
inputStruct.surfaceCase='g'; %Surface type
inputStruct.numPeriods=[1 1 1]; %Number of periods in each direction
inputStruct.levelset=0; %Isosurface level

%% 
% Create stochastic structure

[F,V,C]=triplyPeriodicMinimalSurface(inputStruct);

%%
% Using grouping to keep only largest group
groupOptStruct.outputType='label';
[G,~,groupSize]=tesgroup(F,groupOptStruct); %Group connected faces
[~,indKeep]=max(groupSize); %Index of largest group

%Keep only largest group
F=F(G==indKeep,:); %Trim faces
C=C(G==indKeep,:); %Trim color data 
[F,V]=patchCleanUnused(F,V); %Remove unused nodes

%%
% Get boundary edges 
Eb=patchBoundary(F);

%%
% Visualize surface

cFigure; 
gpatch(F,V,'kw','none');
gpatch(Eb,V,[],'b',1,2);
axisGeom; camlight headlight; 
gdrawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
