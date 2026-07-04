%% quadThick
% Below is a demonstration of the features of the |quadThick| function

%% Syntax
% |[E,VE,Fq1,Fq2]=quadThick(F,V,dirSet,layerThickness,numSteps);|

%% Description
% Use |quadThick| to thicken a quadrilateral mesh to create hexahedral
% elements.  

%%

clear; close all; clc;

%% Examples

%%
% PLOT SETTINGS
fontSize=20;

%% Example: Using |quadThick| to create a hexahedral mesh

% Creating an example polygon
ns=15;
t=linspace(0,pi,ns);
x=cos(t);
y=sin(t);
z=zeros(size(x));
Vc=flipud([x(:) y(:) z(:)]);

%Extruding polygon to a quadrilateral surface
cPar.depth=2; 
cPar.patchType='quad'; 
cPar.dir=0;
cPar.closeLoopOpt=0; 
cPar.numSteps=8;
[F,V]=polyExtrude(Vc,cPar);

%%
% Thickening quadrilaterial surface to hexahedral elements
layerThickness=0.5; 
numSteps=3; 
[E,VE]=quadThick(F,V,1,layerThickness,numSteps);

%Use element2patch to get patch data 
FE=element2patch(E);

%%
% Visualize mesh

cFigure; hold on;
title('Hexahedral mesh');
gpatch(FE,VE,'bw','k',1);
axisGeom;
camlight headlight;
drawnow;

%%
% Create test data set
[F,V]=quadSphere(3,1);

%% 
% Thickening quadrilaterial surface to hexahedral elements
layerThickness=0.4; 
numSteps=3; 
[E,VE,Fq1,Fq2]=quadThick(F,V,-1,layerThickness,numSteps);

%Use element2patch to get patch data 
FE=element2patch(E);

%%
% Visualize mesh

cFigure; hold on;
title('Hexahedral mesh');
gpatch(FE,VE,'bw','k',1);
axisGeom;
camlight headlight;
drawnow;

%%

meshStruct.nodes=VE;
meshStruct.faces=FE;
meshStruct.elements=E;

meshView(meshStruct);

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
