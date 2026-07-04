%% insphere
% Below is a demonstration of the features of the |insphere| function

%%
clear; close all; clc;

%% Syntax
% |[R,Vc]=insphere(E,V);|

%% Description 
% Computes the incentres Vc, and inraddi R of the inspheres for the
% tetrahedral input mesh defined by the element array E and the vertices V.
%% Examples 
% 

%%
% Plot settings

fontSize=20;
faceAlpha1=0.1;

%% 
% % A tetrahedral mesh
boxDim=[2 2 1]; % Box dimenstions
pointSpacing=1; 

[meshStruct]=tetMeshBox(boxDim,pointSpacing);

%%
% Acces output fields
E=meshStruct.elements;
V=meshStruct.nodes;
F=meshStruct.faces;
Fb=meshStruct.facesBoundary;
faceBoundaryMarker=meshStruct.boundaryMarker;

%%
% Compute incentres and inradii for the set of tetrahedra
[R,Vc]=insphere(E,V);

%%
% Visualise inspheres

cFigure; hold on;
title('A tetrahedral mesh showing inspheres','FontSize',fontSize);

hp1=gpatch(F,V,'w','k',faceAlpha1,2);

[Fs,Vs]=geoSphere(2,1);
for q=1:1:size(E,1)
    hp2=gpatch(Fs,Vs.*R(q)+Vc(q,:),'rw','none');
end
legend([hp1,hp2],{'Tetrahedral mesh','inspheres'},'Location','northeastoutside');
axisGeom(gca,fontSize); camlight headlight; 
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
