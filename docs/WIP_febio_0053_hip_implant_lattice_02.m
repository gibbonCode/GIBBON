%% WORK IN PROGRESS

clear; close all; clc;

%%

markerSize1=25;
fontSize=15;

%%

% Load surface geometry
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL');

stlName='hip_implant_iso_merge.stl';
fileName=fullfile(pathName,stlName);
[stlStruct] = import_STL(fileName);
F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};
[F,V]=mergeVertices(F,V);

R=euler2DCM([0 0 0.5*pi]);
V=V*R;

%%
cFigure;
gpatch(F,V,'kw')
axisGeom;
camlight headlight;
drawnow;

%%
[regionA]=tetVolMeanEst(F,V); %Volume for regular tets
stringOpt='-pq1.2AaYQ';
inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=ones(size(F,1),1); %Face boundary markers
inputStruct.regionPoints=[0 0 0]; %region points
inputStruct.regionA=regionA*5;

[meshStruct]=runTetGen(inputStruct); %Run tetGen

%%

E=meshStruct.elements;
Fb=meshStruct.facesBoundary;
V=meshStruct.nodes;

%%
% Visualize input geometry

meshView(meshStruct);

%% Example: Create a dual lattice mesh with outer surface cladding
% See also: |dualClad|

logicLattice_V=V(:,2)<-4.5;
logicLattice_E=all(logicLattice_V(E),2);

FT_other_all=element2patch(E(~logicLattice_E,:),V);
indB=tesBoundary(FT_other_all);
FT_other=FT_other_all(indB,:);

cladOpt=1;
shrinkFactor=0.2;
[FT,VT,CT]=dualLattice(E(logicLattice_E,:),V,shrinkFactor,cladOpt);

%%
% Visualize results

cFigure; hold on;
gpatch(FT_other,V,'w','none',1);
gpatch(FT,VT,CT,'none',1);
% patchNormPlot(FT,VT);
axisGeom;
view(2);
camlight headlight;
axis off;
colormap(viridis(4)); icolorbar;
drawnow;

%%
% Visualize results

cFigure; hold on;
gpatch(FT_other,V,'w','none',1);
gpatch(FT,VT,'w','none',1);
% patchNormPlot(FT,VT);
axisGeom; camlight headlight; view(2);
axis off;
drawnow;

%%

warning('This demo is under construction, no FEBio simulation yet.');

%%
% _*GIBBON footer text*_
%
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
%
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
%
% Copyright (C) 2019  Kevin Mattheus Moerman
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
