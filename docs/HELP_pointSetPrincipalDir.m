%% pointSetPrincipalDir
% Below is a demonstration of the features of the |pointSetPrincipalDir| function

%%
clear; close all; clc;

%% Syntax
% |[V,S,U]=pointSetPrincipalDir(X)|

%% Description

%% Examples

%%
% Plot settings
fontSize=15;
markerSize=25;

%% Using |pointSetPrincipalDir| to determine main directions of a point cloud

%%
% Simulating an ellipsoid with known directions

% Ellipsoid axis stretch factors
ellipStretchTrue=[pi 2 1]

% Create ellipsoid patch data
[F,V,~]=geoSphere(3,1);

V=V.*ellipStretchTrue(ones(size(V,1),1),:);

%Create Euler angles to set directions
E=[-0.25*pi 0.25*pi 0.25*pi]; 
[R_true,~]=euler2DCM(E); %The true directions for X, Y and Z axis
V=(R_true*V')'; %Rotate polyhedron

%%
% This is the true axis system
R_true

%%
% Determine principal directions of the point set (in this case an
% ellipsoidal polyhedron). 

[R_fit]=pointSetPrincipalDir(V)

%% 
% Visualizing results

cFigure; hold on;
gtitle('The polyhedron with true (transparant) and determined (solid) axis directions',fontSize);
gpatch(F,V,0.75*ones(1,3),'none',0.5);
plotV(V,'k.','MarkerSize',markerSize);
hl1=quiverTriad(mean(V,1),R_true,10,[],0.2);
hl2=quiverTriad(mean(V,1),R_fit,7,[],1);
legend([hl1 hl2],{'True directions','Fitted directions'});
camlight('headlight');
axisGeom(gca,fontSize); 
drawnow;

%% Output of the singular value decomposition data data

[R_fit,S,U]=pointSetPrincipalDir(V);

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
