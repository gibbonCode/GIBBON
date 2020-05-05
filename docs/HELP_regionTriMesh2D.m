%% regionTriMesh2D
% Below is a basic demonstration of the features of the |regionTriMesh2D| function.

%% 

clear; close all; clc;

%% CREATING A REGION MESH

%% 
% Creating boundary curves 
V1=[-1 -1; -1 1; 1 1; 1 -1];

%% 
% Meshing the region

% The input variable regionCell is a cell array containing all the boundary
% curves, e.g. for a two curve region 1 we would have something like
% regionSpec{1}={V1,V2} where V1 and V2 are the boundary curves. Multiple
% curves may be given here. The first curve should form the outer boundary
% of the entire region, the curves that follow should define holes inside
% this boundary and the space inside them is therefore not meshed. 

%Defining a region
regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)
plotOn=1; %This turns on/off plotting
pointSpacing=0.1; %Desired point spacing
resampleCurveOpt=1; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,plotOn);
plotV(V1,'b-','LineWidth',2);
axis tight; 
drawnow;

%% CREATING A REGION MESH WITH HOLES

%% 
% Creating boundary curves 

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns+1);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
V1=[x(:) y(:)];

%Boundary 2
[x,y] = pol2cart(t,ones(size(t)));
V2=[x(:) y(:)+4];

%Boundary 3
[x,y] = pol2cart(t,2*ones(size(t)));
V3=[x(:) y(:)-0.5];

%%
% Meshing the region

% The input variable regionCell is a cell array containing all the boundary
% curves, e.g. for a two curve region 1 we would have something like
% regionSpec{1}={V1,V2} where V1 and V2 are the boundary curves. Multiple
% curves may be given here. The first curve should form the outer boundary
% of the entire region, the curves that follow should define holes inside
% this boundary and the space inside them is therefore not meshed. 

%Defining a region
regionCell={V1,V2,V3}; %A region between V1 and V2 (V2 forms a hole inside V1)
plotOn=1; %This turns on/off plotting
pointSpacing=0.5; %Desired point spacing
resampleCurveOpt=1; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,plotOn);

plotV(V1,'b-','LineWidth',2);
plotV(V2,'b-','LineWidth',2);
plotV(V3,'b-','LineWidth',2);
axis tight; 
drawnow;

%% Using input structure instead

%%
% Creating boundary curves 

%Boundary 1
ns=500;
t=linspace(0,2*pi,ns+1);
t=t(1:end-1);
r=5;
a=2;
R=r-(a.*cos(7*(t-pi).^2)-a);
[x,y] = pol2cart(t,R);
V1=[x(:) y(:)];

%Boundary 2
[x,y] = pol2cart(t,(0.75*r)*ones(size(t)));
V2=[x(:) y(:)];

%%
% Meshing the region

%Defining input structure
inputStructure.regionCell={V1,V2};
inputStructure.pointSpacing=0.25; 
inputStructure.resampleCurveOpt=1; 
inputStructure.plotOn=0;

[F,V,boundaryInd]=regionTriMesh2D(inputStructure);

%%

cFigure; hold on;
gpatch(F,V,'r');
plotV(V(boundaryInd,:),'b.','markerSize',15);

axisGeom;
view(2);
drawnow; 

%% Using must points

% Create example boundary curve
V=batman(150);

% Create example interior points
t=linspace(0,2*pi,15)'; t=t(1:end-1);
Vm=[0.4*cos(t) 0.15*sin(t)+0.25];

inputStructure.regionCell={V};
inputStructure.pointSpacing=[]; 
inputStructure.resampleCurveOpt=0; 
inputStructure.plotOn=0;
inputStructure.interiorPoints=Vm;
inputStructure.smoothIterations=250;

[F,V,boundaryInd,interiorInd]=regionTriMesh2D(inputStructure);

%%

cFigure; hold on;
hp(1)=gpatch(F,V,'bw','b',1,2);
hp(2)=plotV(V(boundaryInd,:),'b.','markerSize',25);
% hp(3)=plotV(V(interiorInd,:),'r.','markerSize',35);
axisGeom; 
camlight headlight; 
view(2)
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
