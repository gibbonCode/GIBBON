%% hexMeshBox
% Below is a demonstration of the features of the |hexMeshBox| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=20;
faceAlpha1=1;
faceAlpha2=0.5;

%% CREATING A MESHED BOX
boxDim=[5 6 7];
boxEl=[4 5 6];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E=meshStruct.E;
V=meshStruct.V;
F=meshStruct.F;
Fb=meshStruct.Fb;
faceBoundaryMarker=meshStruct.faceBoundaryMarker;

%%
% Plotting model
cFigure;
subplot(1,2,1);
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha1);
% [hp]=patchNormPlot(Fb,V,1); %Display face normals

axisGeom(gca,fontSize); 
colormap(gjet(6)); icolorbar; 


subplot(1,2,2);
title('Box hexahedral mesh','FontSize',fontSize);
hold on;

gpatch(F,V,0.5*ones(1,3),'k',faceAlpha2);
patchNormPlot(F,V);
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
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
