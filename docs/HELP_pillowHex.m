%% pillowHex
% Below is a demonstration of the features of the |pillowHex| function

%% Syntax
% |[Ep,Vp,Cp]=pillowHex(E,V,C,shrinkFactor);|

%% Description
%

%%
clear; close all; clc;

%%
% Plot settings
fontSize=25;
faceAlpha1=0.25;
edgeWidth=2;
markerSize=35;
cMap=gjet(6);

%% Examples
%

%% Example: Pillowing a hexahedral element

%%
% Creating an example hexahedral element
V=[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1;]; %nodes
E=1:8; %Element

%%

shrinkFactor=0.5; 
[Ep,Vp]=pillowHex(E,V,[],shrinkFactor);

%%
% Visualize results

[F]=element2patch(E);  %Patch data for plotting
[Fp]=element2patch(Ep);  %Patch data for plotting

cFigure;

subplot(1,2,1); hold on;
title('Original hex element','FontSize',fontSize);
gpatch(F,V,'gw','g',0.5,edgeWidth);
% patchNormPlot(F,V,0.25);
plotV(V,'k.','MarkerSize',markerSize);
colormap(cMap);
axisGeom;
axis off;
camlight headlight; 

subplot(1,2,2); hold on;
title('Pillowed hex element','FontSize',fontSize);
gpatch(Fp,Vp,'bw','b',0.5,edgeWidth);
% patchNormPlot(Fp,Vp,0.25);
plotV(Vp,'k.','MarkerSize',markerSize);
colormap(cMap);
axisGeom;
axis off;
camlight headlight; 

drawnow; 

%% Example: Pillowing a set of hexahedral elements

%%
% Creating an example hexahedral element set
boxDim=[5 6 7];
boxEl=[4 5 6];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E=meshStruct.E;
V=meshStruct.V;
C=(1:1:size(E,1))';

%%
% 
shrinkFactor=0.5; 
[Ep,Vp,Cp]=pillowHex(E,V,C,shrinkFactor);

%%
% Visualize results
[F,CF]=element2patch(E,C);  %Patch data for plotting
[Fp,CFp]=element2patch(Ep,Cp);  %Patch data for plotting

cFigure;

subplot(1,2,1); hold on;
title('Original hex element','FontSize',fontSize);
gpatch(F,V,CF,'k',0.5);
% patchNormPlot(F,V,0.25);
plotV(V,'k.','MarkerSize',markerSize);
colormap(cMap);
axisGeom;
axis off;
camlight headlight; 

subplot(1,2,2); hold on;
title('Pillowed hex element','FontSize',fontSize);
gpatch(Fp,Vp,CFp,'k',0.5);
% patchNormPlot(Fp,Vp,0.25);
plotV(Vp,'k.','MarkerSize',markerSize);
colormap(cMap);
axisGeom;
axis off;
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
