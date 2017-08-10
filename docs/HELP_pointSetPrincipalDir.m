%% HELP_pointSetPrincipalDir
% Below is a demonstration of the features of the |pointSetPrincipalDir| function

%%
clear; close all; clc;

%%
% Plot settings

fontSize=11;

%% Using |pointSetPrincipalDir| to determine main directions of a polyhedron

%%
% Simulating an ellipsoid with known directions

% Ellipsoid axis stretch factors
ellipStretchTrue=[pi 2 1]

% Create ellipsoid patch data
[F,V,~]=geoSphere(2,1);

V=V.*ellipStretchTrue(ones(size(V,1),1),:);

%Create Euler angles to set directions
E=[0.25*pi 0.25*pi 0.25*pi]; 
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

cFigure;
title('The polyhedron with true (transparant) and determined (solid) axis directions','FontSize',fontSize);
hold on;

gpatch(F,V,0.75*ones(1,3),'k',1);
quiverTriad(mean(V,1),R_fit,7,[],1);
quiverTriad(mean(V,1),R_true,10,[],0.2);

camlight('headlight');
axisGeom(gca,fontSize); 
drawnow;

%%
% What is clear from the above is that a different system is obtained. This
% is due to the symmetry properties of the ellipsoid. However all vectors
% are colinear with the true vector directions. The output direction matrix
% is ordered in size (as per the singular value decomposition). The vectors
% turned out colinear with R_true due to the fact that the ellipsoid
% directions were biased in a similar sense. However if the order is
% altered the first, second and third axes no longer allign with what was
% viewed here as the true directions. However the singular values can also
% be requested as an output allowing the user to reorder the output
% direction matrix if desired. 

[R_fit,S]=pointSetPrincipalDir(V);
S

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
