%% subTriLocal
% Below is a demonstration of the features of the |subTriLocal| function

%%
clear; close all; clc; 

%% Syntax
% |[outputStruct]=subTriLocal(inputStruct);|

%% Description
% This function refines (using the subTri function) the triangles in F
% defined by indFaces and splits neighbouring triangles to create a
% conforming mesh. This way local refinement can be achieved. 
% The output mesh is stored in Fq and Vq and a face color list Cq can also
% be requested which lists whether a triangle is unaltered (Cq==1), is
% subdevided (Cq==2) or has been split to connect the two regions (Cq==3).
% The optional input f (default is 0) defines the location of the new
% points introduced for the transition elements. Using f>0 (and <1) will
% place these points closer to the coarse mesh nodes. The optional output
% indInitial is a list containing all the original nodes. 

%% Examples

%%
% Plot settings
fontSize=15;
faceAlpha=1;
edgeColor=0.*ones(1,3);
edgeWidth=2;

%% Example: Refining a local region of a mesh (e.g. the top half of a sphere)

%% 
% Building example geometry

%Defining geodesic dome
r=1; %sphere radius
n=2; %Refinements   
[F,V,~]=geoSphere(n,r);

%%
% Define face list for refinement

L=V(:,3)>0.5;
LF=all(L(F),2);
indFaces=find(LF);

inputStruct.F=F; 
inputStruct.V=V; 
inputStruct.indFaces=indFaces; 
[outputStruct]=subTriLocal(inputStruct);
Fq=outputStruct.F; 
Vq=outputStruct.V; 
Cq=outputStruct.faceTypeLabel;

%%
% Plotting results

%Create face color data to visualize selection
C=ones(size(F,1),1);
C(indFaces)=2;

hf=cFigure; 
subplot(1,2,1); hold on; 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

subplot(1,2,2); hold on; 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',Cq,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
colormap(autumn); colorbar; 
view(3); axis tight;  axis equal;  grid on;

%%
% For the above a default f of 0 was assumed. Note the difference when
% instead f=0.25 is used. 
f=0.25; 

inputStruct.F=F; 
inputStruct.V=V; 
inputStruct.indFaces=indFaces; 
inputStruct.f=f; 
[outputStruct]=subTriLocal(inputStruct);
Fq=outputStruct.F; 
Vq=outputStruct.V; 
Cq=outputStruct.faceTypeLabel;

%%
% Plotting results

%Create face color data to visualize selection
C=ones(size(F,1),1);
C(indFaces)=2;

hf=cFigure; 
subplot(1,2,1); hold on; 
title('Initial','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',Cq,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
colormap(autumn); colorbar; 
view(3); axis tight;  axis equal;  grid on;

%% Smoothening the output mesh
% The output mesh can be smoothening normally or using constrained
% smoothing and using the optional indInitial output. The points defined by
% indInitial are plotted on the mesh. Note however that smoothening may
% undo the change induced by the factor f. 

indFaces=1:60;

inputStruct.F=F; 
inputStruct.V=V; 
inputStruct.indFaces=indFaces; 
inputStruct.f=f; 
[outputStruct]=subTriLocal(inputStruct);
Fq=outputStruct.F; 
Vq=outputStruct.V; 
Cq=outputStruct.faceTypeLabel;
indInitial=outputStruct.indInitial;

smoothPar.Alpha=0.1;
smoothPar.Beta=0.5;
smoothPar.Method='HC';
smoothPar.n=100;
smoothPar.RigidConstraints=indInitial;    
[Vq]=tesSmooth(Fq,Vq,[],smoothPar);

%%
% Plotting results

%Create face color data to visualize selection
C=ones(size(F,1),1);
C(indFaces)=2;

hf=cFigure; 
subplot(1,2,1); hold on; 
title('Initial','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',Cq,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
plotV(Vq(indInitial,:),'k.','MarkerSize',50);
camlight headlight;
set(gca,'FontSize',fontSize);
colormap(autumn); colorbar; 
view(3); axis tight;  axis equal;  grid on;

%% Keeping track of face data or color information

%Create example color information
[C]=vertexToFaceMeasure(F,V(:,2));

inputStruct.F=F; 
inputStruct.V=V; 
inputStruct.C=C; 
inputStruct.indFaces=indFaces; 
[outputStruct]=subTriLocal(inputStruct);
Fq=outputStruct.F; 
Vq=outputStruct.V; 
Cq=outputStruct.C;

%%
% Plotting results

hf=cFigure; 
subplot(1,2,1); hold on; 
title('Initial','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','CData',Cq,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(Fq,Vq,0.2);
camlight headlight;
set(gca,'FontSize',fontSize);
colormap(hsv(size(Cq,1))); colorbar; 
view(3); axis tight;  axis equal;  grid on;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>