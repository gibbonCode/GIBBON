%% stochasticMicrostructure
% Below is a demonstration of the features of the |stochasticMicrostructure| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C]=stochasticMicrostructure(M,IND,ptype);|

%% Description
% This function generates Stochastic Bicontinuous Microstructures
%
% Input structure and default values:
%
%   inputStruct.L=1; % characteristic length
%   inputStruct.Ns=80; % number of sampling points
%   inputStruct.Nw=120; % number of waves
%   inputStruct.q0=55; % wave number
%   inputStruct.relD=0.5; % relative density
%   inputStruct.anisotropyFactors=[1 1 1]; %Anisotropy factors
%   inputStruct.isocap=1; %Option to cap the isosurface
%
% Based on: Soyarslan et al. "3D stochastic bicontinuous
% microstructures: Generation, topology and elasticity"
% https://doi.org/10.1016/j.actamat.2018.01.005 
%
% Original author: Sebastien Callens, September 2020

%% Examples

%%
% Plot settings
cMap=parula(250);
faceAlpha1=1;
faceAlpha2=0.5;
edgeColor1='none';
edgeColor2='none';
fontSize=15; 

%% Example 1
inputStruct.L=1; % characteristic length
inputStruct.Ns=75; % number of sampling points
inputStruct.Nw=120; % number of waves
inputStruct.q0=55; % wave number
inputStruct.relD=0.55; % relative density
inputStruct.anisotropyFactors=[1 1 1]; %Anisotropy factors
inputStruct.isocap=1; %Option to cap the isosurface

%% 
% Create stochastic structure

[F,V,C]=stochasticMicrostructure(inputStruct);

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
axisGeom; camlight headlight; 
colormap gjet; icolorbar;
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
