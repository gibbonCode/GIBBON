%% disp2strain
% Below is a demonstration of the features of the |disp2strain| function

%%
clear; close all; clc;

%% Syntax
% |[D_out]=disp2strain(Ux, Uy, Uz, v, strain_type);|

%% Description 
% Computes deformation metrics (including strains) for the meshgrid
% formatted displacement arrays Ux, Uy, and Uz, based on the point spacing
% v and strain_type requested. 

%% Examples 

n=5;
r=linspace(0,1,n); 
v=mean(diff(r))*ones(1,3);
[X1,Y1,Z1]=meshgrid(r);

%Stretches
l1=1.1; l2=0.7; l3=1.5; 

%Scale/stretch coordinates
X2=X1.*l1; Y2=Y1.*l2; Z2=Z1.*l3;

%Get displacements 
Ux=X2-X1; Uy=Y2-Y1; Uz=Z2-Z1;

%%
% Compute deformation data from displacements

% 1 = Biot (linear) strain tensor
% 2 = Hencky (logarithmic/natural) strain tensor
% 3 = Green-Lagrange strain tensor
strain_type=1; 
[D_out]=disp2strain(Ux, Uy, Uz, v, strain_type)

%%

cFigure; 
subplot(1,2,1); hold on; 
plot3(X1(:),Y1(:),Z1(:),'k.','MarkerSize',25); 
axisGeom; 

subplot(1,2,2); hold on; 
plot3(X2(:),Y2(:),Z2(:),'r.','MarkerSize',25); 
axisGeom; 

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
