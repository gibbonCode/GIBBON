%% hexMeshSphere
% Below is a demonstration of the features of the |hexMeshSphere| function

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=1;
edgeColor=0.25*ones(1,3);
edgeWidth=2;

%% Creating a hollow hexahedral mesh sphere 
% Creating a solid hexahedral mesh sphere

%Control settings
cPar.sphereRadius=10;
cPar.coreRadius=5;
cPar.numElementsMantel=5; 
cPar.numElementsCore=8; 
cPar.makeHollow=0;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

%%
% Plotting sphere model

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,~]=element2patch(E(L,:),[],'hex8');

cFigure;
subplot(1,2,1); hold on;
title('The hexahedral mesh sphere','FontSize',fontSize);
gpatch(Fb,V,'r');
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on;
title('Cut-view of the mesh','FontSize',fontSize);
gpatch(Fs,V,'r');
axisGeom(gca,fontSize);
camlight headlight;

drawnow; 

%% Creating a solid hexahedral mesh sphere 
% Creating a solid hexahedral mesh sphere

%Control settings
cPar.sphereRadius=10;
cPar.coreRadius=5;
cPar.numElementsMantel=5; 
cPar.numElementsCore=8; 
cPar.makeHollow=1;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
E=meshStruct.E; %The elements 
V=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces
faceBoundaryMarker=meshStruct.faceBoundaryMarker; %Boundary marker

%%
% Plotting sphere model

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,~]=element2patch(E(L,:),[],'hex8');

cFigure;
subplot(1,2,1); hold on;
title('The hexahedral mesh sphere','FontSize',fontSize);
gpatch(Fb,V,'r');
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on;
title('Cut-view of the mesh','FontSize',fontSize);
gpatch(Fs,V,'r');
axisGeom(gca,fontSize);
camlight headlight;

drawnow; 

%%

cFigure;
title('Boundary color visualization','FontSize',fontSize);
gpatch(Fb,V,faceBoundaryMarker,'none',0.5);
axisGeom(gca,fontSize);
camlight headlight;
colormap(gjet(4)); icolorbar;
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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
