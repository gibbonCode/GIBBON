%% hexMeshHemiSphere
% Below is a demonstration of the features of the |hexMeshHemiSphere| function

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=1;
edgeColor=0.25*ones(1,3);
edgeWidth=2;

%% Creating a heme-sphere hexahedral mesh
% Creating a  heme-sphere hexahedral mesh

%Control settings
optionStruct.sphereRadius=1;
optionStruct.coreRadius=optionStruct.sphereRadius/2;
optionStruct.numElementsMantel=6;
optionStruct.numElementsCore=optionStruct.numElementsMantel*2;
optionStruct.outputStructType=2;
optionStruct.makeHollow=0;
optionStruct.cParSmooth.n=25;

% %Creating sphere
[meshOutput]=hexMeshHemiSphere(optionStruct);

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
E=meshOutput.elements;

%%
% Visualize mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
title('Boundary surfaces','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);
% patchNormPlot(Fb,V);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb,V,'kw','none',0.25);
meshView(meshOutput,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

%% Creating a hollow hexahedral mesh hemisphere 
% Creating a hollow hexahedral mesh hemisphere

%Control settings
optionStruct.sphereRadius=1;
optionStruct.coreRadius=optionStruct.sphereRadius/2;
optionStruct.numElementsMantel=6;
optionStruct.numElementsCore=optionStruct.numElementsMantel*2;
optionStruct.outputStructType=2;
optionStruct.makeHollow=1;
optionStruct.cParSmooth.n=25;

% %Creating sphere
[meshOutput]=hexMeshHemiSphere(optionStruct);

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
E=meshOutput.elements;

%%
% Visualize mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
title('Boundary surfaces','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);
% patchNormPlot(Fb,V);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb,V,'kw','none',0.25);
meshView(meshOutput,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

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
