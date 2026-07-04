%% discQuadMesh
% Below is a demonstration of the features of the |discQuadMesh| function

%% Syntax
% |[varargout]=discQuadMesh(nElements,r,f);|

%% Description 
% 
%% Examples 
% 
clear; close all; clc; 

%% PLOT SETTINGS

fontSize=15;
edgeWidth=2; 

%% Building a quadrilateral circular mesh

ne=4; %Elements in radius
r=2; %Outer radius of disc
f=0.5; %Fraction (with respect to outer radius) where central square appears

%Create the mesh
[Fc,Vc]=discQuadMesh(ne,r,f);

%%
%Visualizing mesh

cFigure; hold on;
title('A circle meshed with quadrilateral elements','FontSize',fontSize);
gpatch(Fc,Vc,'gw','k',1,edgeWidth);
axisGeom(gca,fontSize);
view(2); 
drawnow;

%% Option face label output
% The additional outputs Cc and indEdge provide a color labelling for the
% central square and other elements, and the indices for the outer edge
% respectively. 

[Fc,Vc,Cc,indEdge]=discQuadMesh(ne,r,f);

%%
%Visualizing mesh

cFigure; hold on;
title('A circle meshed with quadrilateral elements','FontSize',fontSize);
gpatch(Fc,Vc,Cc,'k',1,edgeWidth);
axisGeom(gca,fontSize);
view(2); 
drawnow;

%% Constrained smoothing of mesh

Eb=patchBoundary(Fc);

cPar.n=100;
cPar.Method='LAP';
cPar.RigidConstraints=unique(Eb(:));
[Vcs]=patchSmooth(Fc,Vc,[],cPar);

cFigure; 
subplot(1,2,1); hold on;
title('Raw','FontSize',fontSize);
gpatch(Fc,Vc,Cc,'k',1,edgeWidth);
axisGeom(gca,fontSize);
view(2); 

subplot(1,2,2); hold on;
title('Smoothed','FontSize',fontSize);
gpatch(Fc,Vcs,Cc,'k',1,edgeWidth);
axisGeom(gca,fontSize);
view(2); 

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
