%% rhombicDodecahedronMesh
% Below is a demonstration of the features of the |rhombicDodecahedronMesh| function

%%
clear; close all; clc;

%% Syntax
% |[Fc_Q,Fc_T,Ft_Q,Ft_T,Ct_Q,Ct_T,Vt]=rhombicDodecahedronMesh(r,nCopies)|

%% Description
% Creates a rhombic dodecahedron mesh where r sets the radias and nCopies
% (a 1x3 vector) sets the number of copies in the x, y, and z direction.
% The output consists of:
%
% Fc_Q, Fc_T: the quadrilateral and triangular face cell arrays (1 cell
% entry per element). 
%
% Ft_Q, Ft,T: the quadrilateral and triangular face arrays
%
% Ct_Q, Ct,T: color/label data for the face arrays
%
% Vt: the vertex array

%% 
% Plot settings
fontSize=15;
faceAlpha1=1;

%% Creating a mesh of rhombic dodecahedra

r=0.5; %Radii, results in a width of 1
nCopies=[3 3 3]; %Number of offset copies

[E,V,C,F,CF]=rhombicDodecahedronMesh(r,nCopies);

%%
% Plotting results
cFigure; hold on;
title('A mesh of rhombicDodecahedra','FontSize',fontSize);
gpatch(F,V,CF,'k',faceAlpha1);
colormap(gjet); 
axisGeom(gca,fontSize);
camlight('headlight'); 
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
