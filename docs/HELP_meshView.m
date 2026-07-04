%% meshView
% Below is a demonstration of the features of the |meshView| function

%% Syntax
% |[hFig,hp]=meshView(meshStruct,optionStruct);|

%% Description
% This function visualizes volumetric mesh data (e.g. tetrahedral or
% hexahedral elements). The visualization uses the |anim8| function to
% generate a slider controlled cut view of the mesh. The input is a mesh
% structure and an option structure. The default option structure (if
% optionStruct is left empty/incomplete) is:
%  
%   defaultOptionStruct.hFig=[]; %Figure handle (if empty a new figure is
%   created)
%   defaultOptionStruct.numSLiceSteps=25; %Number of slice steps
%   defaultOptionStruct.cMap=[]; %Colormap used (if empty it is based of
%   number of element material types
%   defaultOptionStruct.faceAlpha1=0.2; %Alpha level for boundary surface
%   defaultOptionStruct.faceAlpha2=1; %Alpha level for mesh elements
%   defaultOptionStruct.lightWeightPlot=1; %Option to only plot element
%   outer boundaries to create a more lightweigth plot
%
% The output consists of a figure handle and a graphics handle. 
% If tetgen is used the meshStruct can be the tetgen output structure. 

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.3;
faceAlpha2=1;
cMap=gjet(4); 
patchColor=cMap(1,:);
markerSize=10; 

%% Creating an example mesh

%%
% Surface mesh

testCase=2;
switch testCase
    case 1
        [F,V,~]=geoSphere(2,1); % Building a geodesic dome surface model
    case 2
        [F,V]=stanford_bunny('g'); %Bunny
        V_mean=mean(V,1);
        V=V-V_mean(ones(size(V,1),1),:);
end

%%
% Using |runTetGen| to mesh the geometry

inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=fliplr(F);
inputStruct.Nodes=V;
inputStruct.holePoints=[]; %holes
inputStruct.faceBoundaryMarker=ones(size(F,1),1); %Face boundary markers
inputStruct.regionPoints=getInnerPoint(F,V); %region interior point
inputStruct.regionA=tetVolMeanEst(F,V); %Volume attribute
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput,[]);

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
