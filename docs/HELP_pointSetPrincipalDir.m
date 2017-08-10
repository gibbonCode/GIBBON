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
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
