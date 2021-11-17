%% DEMO_mesh_bifurcation_angle_control
% This demo shows the use of the |splitCurveSetMesh| function to
% parameterise a bifurcation in terms of vessel directions and diameters.

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=15;
lineWidth=3;
markerSize1=25;

%% Control parameters

r1=3;
r=(0.5*r1^3)^(1/3);
r2=r*0.9;
r3=r*0.75;
pointSpacingMain=0.5;
V1_origin=[0 0 0]; %Origin of first circle
bifurcationAngleDeg2=35;
bifurcationAngleDeg3=-45;
bifurcationDistance2=3;
bifurcationDistance3=3;

nSmooth=25; %Number of Laplacian/HC smoothing steps
splitMethod='ortho'; %'nearMid'; %Saddle placement 
w=1; %1=saddle arcs upward to max height, 0 means saddle is in plane of first circle

%% Derived metrics

nz=[0 0 1]; %Normal direction for z-axis
R2=euler2DCM([0 (bifurcationAngleDeg2/180)*pi 0]); %Rotation matrix for first direction
R3=euler2DCM([0 (bifurcationAngleDeg3/180)*pi 0]); %Rotation matrix for second direction

n2=nz*R2; %Direction vector for first branch
n3=nz*R3; %Direction vector for second branch

V2_origin=n2.*bifurcationDistance2; %Origin of second circle
V3_origin=n3.*bifurcationDistance3; %Origin of third circle

%Number of points to use allong circle 1
np=ceil((2*pi*r1)./pointSpacingMain);
np=np+~iseven(np); %Forcing this even creates symmetric saddle position

%% Create curves

% Circle 1
t=linspace(0,2*pi,np+1)'; t=t(1:end-1);
x=r1.*sin(t(:));
y=r1.*cos(t(:));
z=zeros(size(t));
V1=[x y z];
V1=V1+V1_origin;

% Circle 2
x=r2.*sin(t);
y=r2.*cos(t);
z=zeros(size(x));
V2=[x y z]*R2;
V2=V2+V2_origin+V1_origin;

% Circle 3
x=r3.*sin(t);
y=r3.*cos(t);
z=zeros(size(t));
V3=[x y z]*R3;
V3=V3+V3_origin+V1_origin;

%% Meshing bifurcation

numStepsBranch=ceil((bifurcationDistance2+bifurcationDistance3)/2./pointSpacingMain); 
V_cell={V1,V2,V3};
patchType='quad';
smoothPar.Method='HC';
smoothPar.n=nSmooth;

[F,V,curveIndices,faceMarker]=splitCurveSetMesh(V_cell,numStepsBranch,patchType,smoothPar,splitMethod,w);

%% Visualization

cFigure; 
subplot(1,2,1); hold on;
title('Input contours and output mesh','FontSize',fontSize);
gpatch(F,V,faceMarker);
plotV(V(curveIndices{1},:),'r.-','MarkerSize',markerSize1,'LineWidth',lineWidth);
plotV(V(curveIndices{2},:),'b.-','MarkerSize',markerSize1,'LineWidth',lineWidth);
plotV(V(curveIndices{3},:),'y.-','MarkerSize',markerSize1,'LineWidth',lineWidth);
quiverVec(V1_origin,n2,bifurcationDistance2,'b')
quiverVec(V1_origin,n3,bifurcationDistance3,'y')
axisGeom(gca,fontSize);
colormap parula; icolorbar;
camlight headlight; 

subplot(1,2,2); hold on;
title('Output mesh','FontSize',fontSize);
gpatch(F,V,'w','k');
axisGeom(gca,fontSize);
camlight headlight; 

drawnow;

% %% Smooth Catmull-Clark subdevision 
% 
% fixBoundaryOpt=1;
% nIter=2;
% 
% [Fs,Vs]=subQuadCatmullClark(F,V,nIter,fixBoundaryOpt);
% 
% cFigure; hold on;
% title('Input contours and connected mesh','FontSize',fontSize);
% gpatch(Fs,Vs,'w','k');
% axisGeom(gca,fontSize);
% camlight headlight; 
% drawnow;

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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
