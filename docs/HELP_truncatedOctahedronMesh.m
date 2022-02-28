%% truncatedOctahedronMesh
% Below is a demonstration of the features of the |truncatedOctahedronMesh| function

%%
clear; close all; clc;

%% Syntax
% |[F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(r,nCopies)|

%% Description
% Creates a truncated octahedron mesh where r sets the radias and nCopies
% (a 1x3 vector) sets the number of copies in the x, y, and z direction.
% The output consists of:
%
% F1c, F2c: the hexgonal and quadrilateral face cell arrays (1 cell entry
% per element).  
%
% F1, F2: the hexgonal and quadrilateral face arrays
%
% C1, C2: color/label data for the face arrays
%
% VT: the vertex array

%% 
% Plot settings
fontSize=15;
faceAlpha1=1;

%%

r=sqrt(5)/4; %Radii, results in a width of 1

n=3; %Desired number of copies in each direction 

%The actual input 
nCopies=[n n n+ceil((n+1)/2)+1]; %Number of offset copies


%%
[F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(r,nCopies);

%%
% Plotting results
cFigure; hold on;
gtitle('A mesh of truncated octahedra');
gpatch(F1,VT,C1,'k',faceAlpha1);
gpatch(F2,VT,C2,'k',faceAlpha1);
colormap(gjet);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow; 

%%
% Plotting results
cFigure; hold on;
gtitle('A mesh of truncated octahedra');
gpatch(F1,VT,'none','k',faceAlpha1,3);
gpatch(F2,VT,'none','k',faceAlpha1,3);
colormap(gjet);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow; 

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
