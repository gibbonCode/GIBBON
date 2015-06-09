%% regionTriMesh3D
% Below is a basic demonstration of the features of the |regionTriMesh3D| function.

%% Syntax
% |[F,V]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);|

%% Description
% This function meshes a 3D region. It is conceptually similar to
% regionTrimMesh2D. The 3D regions are rigidly transformed to a nearly 2D
% form and |regionTrimMesh2D| is then used for 2D meshing. Then
% interpolation is to get the coordinates in the 3rd dimension. 

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
markerSize1=45;
lineWidth1=4;
faceAlpha=0.5;

%% Example: Meshing a 3D region

%%
% Creating example boundary curves and input region specification

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
z=1/10*x.^2;
V1=[x(:) y(:) z(:)];

%Boundary 2
ns=100;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
[x,y] = pol2cart(t,ones(size(t)));
z=zeros(size(x));
V2=[x(:) y(:)+4 z(:)];

%Boundary 3
ns=75;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
[x,y] = pol2cart(t,2*ones(size(t)));
z=zeros(size(x));
V3=[x(:) y(:)-0.5 z(:)];

%Create Euler angles to set directions
E=[0.25*pi -0.25*pi 0];
[R,~]=euler2DCM(E); %The true directions for X, Y and Z axis

V1=(R*V1')'; %Rotate polygon
V2=(R*V2')'; %Rotate polygon
V3=(R*V3')'; %Rotate polygon

regionCell={V1,V2,V3}; %A region between V1 and V2 (V2 forms a hole inside V1)

%% 
% Meshing the region (See also |regionTriMesh2D|)

%Defining a region and control parameters (See also |regionTriMesh2D|)
pointSpacing=0.5; %Desired point spacing
resampleCurveOpt=1; 
interpMethod='linear'; %or 'natural'
[F,V]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);

%%
% Plotting meshed model
hf1=cFigure;
title('The meshed model','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','g');

plotV(V1,'b-','LineWidth',2);
plotV(V2,'b-','LineWidth',2);
plotV(V3,'b-','LineWidth',2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
camlight headlight;
drawnow;

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>