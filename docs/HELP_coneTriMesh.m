%% subtri
% Below is a demonstration of the features of the |subtri| function

%% Syntax
% |[F,V]=coneTriMesh(coneRadius,coneHeight,pointSpacing,nTopMin);

%% Description
% The |coneTriMesh| function creates a triangulated mesh for a cone with a
% base radius coneRadius and a height coneHeight. The desired node spacing
% is set by pointSpacing. The parameter nTopMin sets the number of repeated
% "pie-segments" used to build the cone, and hence also sets the minimum
% number of points the tip node is connected to. 

%% Examples

clear; close all; clc;

%%

coneRadius=2; 
coneHeight=4; 
pointSpacing=0.25;
nTopMin=6;
closeBaseOpt=1;

[F,V,C]=coneTriMesh(coneRadius,coneHeight,pointSpacing,nTopMin,closeBaseOpt);

%%

actualMeanPointSpacing = mean(patchEdgeLengths(F,V))


%%
cFigure; hold on; 
title('A triangulated cone')
hp=gpatch(F,V,C,'w',1,1);
% patchNormPlot(F,V);
% gedge(patchBoundary(F),V,'r',5)
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
