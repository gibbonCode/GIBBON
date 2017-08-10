%% gaborCoil
% Below is a demonstration of the features of the |gaborCoil| function

%%
clear; close all; clc;

%% Syntax
% |[V_coil_rep]=gaborCoil(varargin);|

%% Description 
% UNDOCUMENTED 

%% Examples 

%% Creating a single Gabor coil
% 

V=[0 0 0; 1 0 0];
E=[1 2];
numSteps=250;
numTwist=4; 
coilAmplitude=0;
f=3;
Vg=gaborCoil(V,E,numSteps,numTwist,coilAmplitude,f);

%%
% Visualizing coil curve
cFigure; 
hold on; 
plotV(V,'r.','markerSize',50);
plotV(Vg,'k-','LineWidth',5);
axisGeom; 
drawnow; 

%% Creating Gabor coil on all edges in mesh

%%
% Creating example patch data
[F,V]=geoSphere(1,1);
% [F,V]=stanford_bunny;

%%
% Get patch edges
E=patchEdges(F,1);

%%

numSteps=100;
numTwist=8; 
coilAmplitude=0.05; 
f=3;

%%
cFigure; hold on;
gpatch(F,V,'g',0.5*ones(1,3),0.5);
plotV(V,'r.','MarkerSize',25);
axisGeom;
camlight headlight; 

%%

[V_coil_rep]=gaborCoil(V,E,numSteps,numTwist,coilAmplitude,f);

for q=1:1:size(V_coil_rep,3)
    plotV(V_coil_rep(:,:,q),'k-','lineWidth',5);               
end

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
