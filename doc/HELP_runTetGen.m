%% runTetGen
% Below is a demonstration of the features of the |runTetGen| function
%
%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
patchColor=[1 0.5 0];
markerSize=15; 
cMap=gjet(250); 

%%
% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(fileparts(filePath)),'data','temp');
modelName=fullfile(savePath,'tetgenmodel');

%% MESHING A SINGLE REGION MODEL

%%
% Building a geodesic dome surface model
[F,V,~]=geoSphere(2,1);

%%
% Plotting model
hf=cFigure;
title('Surface model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor',patchColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

%%
% DEFINE FACE BOUNDARY MARKERS
faceBoundaryMarker=ones(size(F,1),1);

%%
% Define region points
V_regions=[0 0 0];

%%
% Define hole points
V_holes=[];

%% 
% Regional mesh volume parameter
[regionA]=tetVolMeanEst(F,V); %Volume for regular tets

%% 
% CREATING THE INPUT STRUCTURE
stringOpt='-pq1.2AaY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% 
% PLOTTING MODEL 

%Selecting half of the model to see interior
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L),'tet4');

hf1=cFigure;
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Cut view of Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor,'Marker','.','MarkerEdgeColor','k','MarkerSize',markerSize);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% VIEWING THE MODEL IN TETVIEW
% TetView is an external (non-MATLAB) program for viewing TetGen meshes.
% The runTetGen function also copies the mesh output files to the
% tetView directory. TetView is usually found here:
% ...\gibbon\trunk\lib_ext\tetGen
% You can run TetView seperately or use the following to view the model in TetView:
[runStatus,runCmdHist]=runTetView(meshOutput.loadNameStruct.loadName_ele);

%% 
% Here is an example screeshot for viewing models in tetView:
% 
% <<tetView_screenshot.jpg>>
%
%% MESHING A MULTI-REGION MODEL

%%
% Simulating a multiregion mesh
[F1,V1]=parasaurolophus; %A dino
[F2,V2,~]=geoSphere(2,0.4); %An internal region
V2(:,1)=2*V2(:,1);
V_centre=[0.75 0 0.25];
V2=V2+V_centre(ones(size(V2,1),1),:);

%Joining surface sets
F=[F1;F2+size(V1,1)];
V=[V1;V2];
C=[ones(size(F1,1),1);2*ones(size(F2,1),1)]; %Surface marker colors

%%
% Plotting model
hf=cFigure;
title('Multi-region surface model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
colormap(cMap); 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

%%
% DEFINE FACE BOUNDARY MARKERS
faceBoundaryMarker=C;

%%
% Define region points
V_regions=[0.5 0 0;V_centre];

%%
% Define hole points
V_holes=[];

%% 
% Regional mesh parameters
[A]=tetVolMeanEst(F,V);
regionA=[A A*3];

%% 
% CREATING THE INPUT STRUCTURE

stringOpt='-pq1.2AaYQ';
inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;


%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% 
% PLOTTING MODEL 

%Selecting half of the model to see interior
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L),'tet4');

hf1=cFigure;
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Cut view of Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor,'Marker','.','MarkerEdgeColor','k','MarkerSize',markerSize);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% MESHING A MULTI-REGION MODEL CONTAINING HOLES 

%%
% Simulating a multiregion mesh

%A dino
[F1,V1]=parasaurolophus; 

%An internal region
[F2,V2,~]=geoSphere(2,0.5); 
V2(:,1)=2*V2(:,1);

%An internal hole
[F3,V3,~]=geoSphere(3,0.35); 

%Centering internal structures
V_centre=[0.75 0 0.25];
V2=V2+V_centre(ones(size(V2,1),1),:);
V3=V3+V_centre(ones(size(V3,1),1),:);

%Joining surface sets
F=[F1;F2+size(V1,1);F3+size(V1,1)+size(V2,1)];
V=[V1;V2;V3];
C=[ones(size(F1,1),1);2*ones(size(F2,1),1);3*ones(size(F3,1),1)]; %Surface marker colors

%%
% Plotting model
hf=cFigure;
title('Multi-region surface model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

hp=patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha1,'edgeColor','none');
camlight headlight;
colormap(cMap); 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

%%
% DEFINE FACE BOUNDARY MARKERS
faceBoundaryMarker=C;

%%
% Define region points

%Find an interior point in region 1
logicRegion=ismember(faceBoundaryMarker,[1 2]);
[V_in_1]=getInnerPoint(F(logicRegion,:),V);
plotV(V_in_1,'r.','MarkerSize',25);

%Find an interior point in region 1
logicRegion=ismember(faceBoundaryMarker,[2 3]);
[V_in_2]=getInnerPoint(F(logicRegion,:),V);
plotV(V_in_2,'r.','MarkerSize',25);

V_regions=[V_in_1; V_in_2];

%%
% Define hole points
V_holes=[V_centre];

%% 
% Regional mesh parameters
[A]=tetVolMeanEst(F,V);
regionA=[A A/2];

%% 
% CREATING THE INPUT STRUCTURE 

stringOpt='-pq1.2AaQY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% 
% PLOTTING MODEL 

%Selecting half of the model to see interior
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L),'tet4');

hf1=cFigure;
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Cut view of Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor,'Marker','.','MarkerEdgeColor','k','MarkerSize',markerSize);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% Meshing from a quadrilateral input surface
% Build a quadrilateral surface

boxDim=[4 5 6];
boxEl=[4 5 6];
[Fq,Vq,faceBoundaryMarker_q]=quadBox(boxDim,boxEl);

%%
% CREATING THE INPUT STRUCTURE

[regionA]=tetVolMeanEst(Fq,Vq); %Volume for regular tets

stringOpt='-pq1.2AaYQ';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=Fq;
inputStruct.Nodes=Vq;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=faceBoundaryMarker_q; %Face boundary markers
inputStruct.regionPoints=[0 0 0]; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% 
% PLOTTING MODEL 

%Selecting half of the model to see interior
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L),'tet4');

hf1=cFigure;
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Cut view of Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor,'Marker','.','MarkerEdgeColor','k','MarkerSize',markerSize);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% Meshing 10-node (i.e. quadratic) tetrahedral elements

% Building spherical surface models
r1=1;
r2=2*(r1/3);
r3=r1/3;
[F_sphere,V_sphere,~]=geoSphere(2,r1);
[F_core,V_core,~]=geoSphere(1,r2);
[F_core2,V_core2,~]=geoSphere(1,r3);

% Merging node sets
V1=[V_sphere;V_core;V_core2];

F1=[F_sphere;F_core+size(V_sphere,1);F_core2+size(V_sphere,1)+size(V_core,1)];

% Face boundary markers
faceBoundaryMarker=[ones(size(F_sphere,1),1); 2*ones(size(F_core,1),1);  3*ones(size(F_core2,1),1)];

%%
% Define region point
V_regions=[0 0 (r1+r2)/2;0 0 (r2+r3)/2;];

% Define hole points
V_holes=[0 0 0];

% Regional mesh parameter (desired volume)
[A]=tetVolMeanEst(F1,V1);
regionA=A*ones(1,2); 

stringOpt='-pq1.2AaYQ';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F1;
inputStruct.Nodes=V1;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%%
% Setting desired element type to tet10 (default is tet4)
inputStruct.tetType='tet10';

%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% 
% PLOTTING MODEL 

%Selecting half of the model to see interior
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L),'tet10');

hf1=cFigure;
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Cut view of Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor,'Marker','.','MarkerEdgeColor','k','MarkerSize',markerSize);
view(3); axis tight;  axis equal;  grid on;
colormap(cMap); 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
