%% edgeLengths
% Below is a demonstration of the features of the |edgeLengths| function

%%
clear; close all; clc;

%% Syntax
% |[D]=edgeLengths(E,V);|

%% Description 
% This function computes the edge lengths D for the edges defined by the
% input edges array E and the vertices V. 
%
% See also: |patchEdgeLengths|

%% Examples 
% 

%%
% Create example data 

%Create an ellipsoid from a sphere mesh (so there are various edge lengths)
[F,V]=geoSphere(1,1); %Geodesic sphere mesh
V(:,1)=2*V(:,1); %Stretch into ellipsoid

%Get edge array 
E=patchEdges(F,V);

%%
% Compute edge lenths
D=edgeLengths(E,V);

%%
% VisualiZing edge length data

cFigure; hold on; 
title('Edge lengths');
gpatch(F,V,'w','none');
gedge(E,V,D,2)
axisGeom; camlight headligth; 
colormap gjet; colorbar; 
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
