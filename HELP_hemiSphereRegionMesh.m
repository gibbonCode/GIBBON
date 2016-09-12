%% hemiSphereRegionMesh
% Below is a demonstration of the features of the |ind2patch| function

%%
close all; clc; clear;

%% Syntax
% |[F,V,regionIndSub]=hemiSphereRegionMesh(hemiSphereStruct);|

%% Description
% The |hemiSphereRegionMesh| function creates the faces (F), vertices (or
% nodes, V) and the region indices (regeionIndSub) for a hemi-sphere
% according to the input structure hemiSphereStruct. The latter defines the
% sphere radius, the number of refinement steps for the regions and the
% number of refinement steps for the mesh. For more information on the
% refinement see the |geoSphere| and |subTri| functions and associated demo
% files. 
% A complete sphere is first represented as an icosahedron which is then
% refined (subtriangulated) hemiSphereStruct.nRefineRegions times (whereby
% for each iteration each triangle is subdevided into 4 triangles). This
% initial subdevision defines the element regions. The next refinement step
% defines the number of triangles for each region. The field
% hemiSphereStruct.nRefineMesh defines how many times each mesh region is
% iteratively subtriangulated. 

%% Examples

%%
% Plot settings
fontSize=25;
faceAlpha=1;
lineWidth=1;
markerSize=5;

%%
close all; clc; clear;

%% Example:  Creating a hemisphere mesh using the |hemiSphereRegionMesh| function
% Defining hemi-sphere parameters
hemiSphereStruct.sphereRadius=1; %Sphere radius
hemiSphereStruct.nRefineRegions=2; %Number of refinement steps for regions
hemiSphereStruct.nRefineMesh=2; %Number of refinement steps for mesh

% Get hemi-sphere mesh
[F,V,regionIndSub]=hemiSphereRegionMesh(hemiSphereStruct);

%% 
% Plotting results
%Creating a random color for the each mesh region
cmap=hsv(max(regionIndSub(:)));
cmap=cmap(randperm(size(cmap,1)),:); %scramble colors

hf=cFigure;
hold on; view(3); 
title('Half dome showing regions with subtriangulated mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','flat','EdgeColor','k','CData',regionIndSub,'FaceAlpha',faceAlpha,'LineWidth',lineWidth,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',markerSize);

colormap(cmap); %colorbar; 
axis tight;  axis equal;  grid on;
set(gca,'FontSize',fontSize);
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
