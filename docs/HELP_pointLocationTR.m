%% pointLocationTR
% Below is a demonstration of the features of the |pointLocationTR| function

%%
clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceColor1='g';
faceColor2='r';
faceAlpha1=0.3;
faceAlpha2=1;
edgeColor=0.4*ones(1,3);
edgeWidth=2;
markerSize=2;
cMap=jet(250);

%% FINDING POINTS INSIDE A TESSELATION

%%
% Create a test tesselation
r=1; %Radius of tetrahedron circumsphere
[V,F]=platonic_solid(1,r);

E=[1 2 3 4];

n=3; 
E=E; V=V; 
for q=1:1:n
    [E,V]=subTet(E,V,1);
end
C=(1:size(E,1))';

%%
% Define query points
n=50; 
[X,Y,Z]=meshgrid(linspace(-r,r,n));

QP=[X(:) Y(:) Z(:)];

%%
% Plotting the test case
[F,CF]=element2patch(E,C);    

cFigure;
title(['A ',num2str(size(E,1)),' element tesselation with ',num2str(size(QP,1)),' query points'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);
plotV(QP,'b.','MarkerSize',markerSize);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
view([-50,12])
set(gca,'FontSize',fontSize);
drawnow;

%%
% Convert tesselation to triangulation class
TR = triangulation(E,V);

%% 
% Use |pointLocationTR| to test points
[ti,bc]=pointLocationTR(TR,QP,1,1,1);

%%
% Plot results and color points to the element number they are contained in
cFigure;
title('Interior query points colored to element index','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

logicPlot=~isnan(ti); %Point selection logic
% plotV(QP(logicPlot,:),'b.','MarkerSize',markerSize);
scatter3(QP(logicPlot,1),QP(logicPlot,2),QP(logicPlot,3),75,ti(logicPlot),'fill','markerEdgeColor','k');

colormap(cMap); colorbar;
view(3); axis equal; axis tight; axis off; 
view([-50,12])
set(gca,'FontSize',fontSize);
drawnow;

%% EXAMPLE USING BARY CENTRIC COORDINATE OUTPUT TO MAP POINTS IN A DEFORMED MODEL
% Changing shape
V=TR.Points; 
V(:,1)=V(:,1)+0.5.*sin(pi*V(:,3)); 
V(:,2)=V(:,2)+0.5.*cos(pi*V(:,3)); 

% Fixing triangulation
TR = triangulation(E,V);
QPm = barycentricToCartesian(TR,ti(logicPlot),bc(logicPlot,:));

%%
% Plotting deformed shape
cFigure;
title('Mapping of points in deformed object','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

scatter3(QPm(:,1),QPm(:,2),QPm(:,3),75,ti(logicPlot),'fill','markerEdgeColor','k');

colormap(cMap); colorbar;
view(3); axis equal; axis tight; axis off; 
view([-50,12])
set(gca,'FontSize',fontSize);
drawnow;

%% EXAMPLE USE FOR FINDING VOXELS INSIDE TETRHEDRAL ELEMENTS
%% 
% Simulate image data
load mri;
M=double(squeeze(D)); %example image data set
clear D;
M=M(1:2:end,1:2:end,:);

voxelSize=2./[1,1,.4]; %Voxel size
voxelSize(1)=voxelSize(1)*2;
voxelSize(2)=voxelSize(2)*2;
siz=size(M); %Image size
FOV=siz.*voxelSize; %Field of view size

%% 
% Get patch data
T_low=min(M(:))+((max(M(:))-min(M(:)))/10); %Threshold example
logicVoxels=(M>T_low);
[Fv,Vv,Cv]=ind2patch(logicVoxels,M,'vb'); 
Cv=Cv./max(Cv(:));
[Vv(:,1),Vv(:,2),Vv(:,3)]=im2cart(Vv(:,2),Vv(:,1),Vv(:,3),voxelSize); 

%%
% Simulate a tesselated model

Vv_mean=mean(Vv,1);
[V,F]=platonic_solid(1,max(FOV)/2);
V=V+Vv_mean(ones(1,size(V,1)),:);
E=[1 2 3 4];

n=4; 
E=E; V=V; 
for q=1:1:n
    [E,V]=subTet(E,V,1);
end
C=[1:size(E,1)]';

%%
[F,CF]=element2patch(E,C);    

cFigure;
title('Image data and a tetrahedral tesselation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);
hp1= patch('Faces',Fv,'Vertices',Vv,'FaceColor','flat','FaceVertexCData',Cv(:,ones(1,3)),'EdgeColor','k','FaceAlpha',faceAlpha1);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
drawnow;

%%
% Setting up query points, in this case voxel coordinates
[Jv,Iv,Kv]=meshgrid(1:1:size(M,2),1:1:size(M,1),1:1:size(M,3)); %Voxel image coordinates
[Xv,Yv,Zv]=im2cart(Iv,Jv,Kv,voxelSize); %Voxel cartesian coordinates
Vq=[Xv(:) Yv(:) Zv(:)];

%% 
% Find what elements the voxels are contained in
TR = triangulation(E,V);
[TI,BC]=pointLocationTR(TR,Vq,1,1,1); 

% Reshaping the output TI leads in effect to a labeled image whereby each
% voxels containes an index (or group/label number) for an element
ML=reshape(TI',size(M));

%%
% Plotting the result
cFigure;
title('Image data found inside the tetrahedral tesselation','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','g','EdgeColor','none','FaceAlpha',0.1);

logicVoxels=~isnan(ML);
[Fm,Vm,Cm]=ind2patch(logicVoxels,M,'vb');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha2);

colormap gray; colorbar; 
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% EXAMPLE USE FOR INTERPOLATION OF IMAGE DATA ONTO VOLUME ELEMENTS
% Hence deriving image based intensity values becomes a matter of averaging
% across each element index label. If the voxels are sufficiently small
% with respect to the elements then this method actually partially takes
% into account partial voluming effects. i.e. the intensity of an element
% becomes the mean of the voxels contained within it which in some cases is
% more accurate then an interpolation which does not take into account the
% volumetric nature of the mesh (e.g. methods based on
% sampling/interpolating the intensities at the centre point or nodal
% points of the mesh only instead). 
% 
[E_color]=imlabelMean(M,ML); %The image based element colors

%%
% Plotting results

cFigure;
title('The image data and the tesselation colored towards the image data','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

[F,C]=element2patch(E,E_color); 

L=~isnan(C);
patch('Faces',F(L,:),'Vertices',V,'FaceColor','flat','CData',C(L),'EdgeColor','k','FaceAlpha',faceAlpha2);

L_plot=false(siz);
L_plot(:,:,round(size(M,3)/2))=1;
L_plot=L_plot;
[Fm,Vm,Cm]=ind2patch(L_plot,double(M),'sk');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha2);

L_plot=false(siz);L_plot(round(size(M,1)/2),:,:)=1;
L_plot=L_plot;
[Fm,Vm,Cm]=ind2patch(L_plot,M,'si');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha2);

L_plot=false(siz);L_plot(:,round(size(M,2)/2),:)=1;
L_plot=L_plot;
[Fm,Vm,Cm]=ind2patch(L_plot,M,'sj');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha2);

colormap gray; colorbar; 
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

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
