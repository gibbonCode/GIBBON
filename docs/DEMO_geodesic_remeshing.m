%% DEMO_geodesic_remeshing
% Below is a demonstration for:
% 
% * The use of remeshTriSurfDistMap (and subTri and smoothening) for surface remeshing

%%

clear; close all; clc;

%%
%Plot settings
cMapDist=flipud(igviridis(250));
[cMapIndices,scrambleIndices]=scramble(viridis(250),1); %Colormap

faceAlpha1=1;
faceAlpha2=0.65;
fontSize=25; 
markerSize=50;
lineWidth=4;
scatterSize=65;

%% Example 1: Resampling a triangulated surface model

%% Example: Using |meshDistMarch| for geodesic point sampling

%%
% Get example triangulated mesh data
[F,V]=graphicsModels(7);

%% 
% Visualize input mesh

cFigure; hold on; 
title('Input mesh')
gpatch(F,V,'gw'); 
axisGeom(gca,fontSize); 
camlight headlight;
drawnow; 

%%
% Compute distances on mesh

numSeeds=120;

%Option set
[~,indStart]=min(V(:,1)); %Index of the start point
optionStruct.toleranceLevel=0; %Tolerance for convergence
optionStruct.waitBarOn=1; %Turn on/off waitbar

%Compute distances on mesh description
[Fn,Vn,seedIndex,indSeeds,d]=remeshTriSurfDistMap(F,V,numSeeds,indStart,optionStruct); %distance based marching
[~,~,ind2]=unique(seedIndex);

%%
% Visualization

cFigure;
subplot(1,3,1); hold on;
title('Distances on a triangulated surface model','fontSize',fontSize);
hp(1)=gpatch(F,V,d,'none',1); hp(1).FaceColor='Interp';
hp(2)=plotV(V(indSeeds,:),'k.','MarkerSize',markerSize);
legend(hp,{'Mesh distances','Seed point(s)'},'Location','SouthOutSide');
axisGeom;
camlight headlight;
colormap(gca,cMapDist); colorbar;

subplot(1,3,2); hold on;
title('Seed indices','fontSize',fontSize);
hp(1)=gpatch(F,V,ind2,'none',1); 
hp(2)=plotV(V(indSeeds,:),'k.','MarkerSize',markerSize);
legend(hp,{'Mesh seed indices','Seed point(s)'},'Location','SouthOutSide');
axisGeom;
camlight headlight;
colormap(gca,cMapIndices); 

subplot(1,3,3); hold on;
title('Resampled surface model','fontSize',fontSize);
plotV(V(indSeeds,:),'k.','MarkerSize',50);
hp(1)=gpatch(F,V,'kw','none',0.5);
hp(2)=gpatch(Fn,Vn,'gw','k',1,2);
legend(hp,{'Original mesh','Resampled mesh'},'Location','SouthOutSide');

axisGeom;
camlight headlight;
drawnow;

%% Example 2: Boundary preserving resampling

%%
% Control parameters
pointSpacing=1; %Desired point spacing
nRefineOriginal=1; %Set n>0 to refine the input mesh through sub-triangulation before resampling
nRefineOutput=0; %Number of output refinement steps

%% Creating example surface data

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
z=1/10*x.^2;
V1=[x(:) y(:) z(:)];

%Create Euler angles to set directions
E=[0.25*pi -0.25*pi 0];
[R,~]=euler2DCM(E); %The true directions for X, Y and Z axis

V1=(R*V1')'; %Rotate polygon

regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)

% Meshing the region (See also |regionTriMesh2D|)
[F,V]=regionTriMesh3D(regionCell,0.15,1,'natural');

%% Refine input mesh before resampling
% Refining the input mesh is required if the intended mesh density exceeds
% that of the desired output mesh density. For each iteration the triangle
% edges are split in half while the triangle faces are split into 4. The
% mesh becomes very dense, very quickly so do not over do this. Subtri
% works through triangle splitting and is linear in the sense that it
% leaves input points unaltered but adds intermediate points on all input
% edges. 

if nRefineOriginal>0
    for q=1:1:nRefineOriginal
        [F,V]=subtri(F,V,1); %Refine input mesh through sub-triangulation
    end
end

%% Estimate number of points required given point spacing

numPointsInput=size(V,1); %Number of points in the original data
[A]=patch_area(F,V); %Areas of current faces
totalArea=sum(A(:)); %Total area
l=sqrt(totalArea); %Width or length of square with same size
np=round((l./pointSpacing).^2); %Point spacing for mesh in virtual square

%% Get indices of boundary points and get boundary curve

%Boundary edges
Eb=patchBoundary(F,V);

%Reorder edge list to obtain indices describing continuous curve (assuming
%a single boundary, need to group first and do reordering on each group if
%multiple boundaries exist) 
[indList]=edgeListToCurve(Eb);
indList=indList(1:end-1);

%% Get indices of "must points" close to evenly spaced on boundary curve
 
D=pathLength(V(indList,:)); %The cummulative curve length
boundaryLength=max(D); %The total curve length

nb=round(boundaryLength./pointSpacing); %Number of points to keep on boundary
[Vb] = evenlySampleCurve(V(indList,:),nb,'pchip',1); 

[~,minIND]=minDist(Vb,V(indList,:));
indListSelect=indList(minIND); %List of points to keep

%% Visualize boundary points to keep on input mesh

cFigure; hold on; 
title('Input mesh and boundary points to keep')
gpatch(F,V,'gw'); 

plotV(V(indListSelect,:),'r.-','MarkerSize',25,'LineWidth',lineWidth);
plotV(Vb,'b.-','MarkerSize',markerSize,'LineWidth',lineWidth);

axisGeom(gca,fontSize); 
camlight headlight;
drawnow; 

%% Resample input surface geodesically
% Geodesic re-sampling works by taking a point on the surface and
% calculating the distance to all other points by marching/propagating
% across the mesh. The point on the mesh furthest away from the input point
% (or points) is then added to the list. Then new distances are computed
% and the furthest point is again added to this list. Therefore using this
% iterative process equally spaced point sets can be obtained. The mesh
% works by sampling the input mesh to a coarser homogeneous mesh. If the
% input mesh is too coarse refine it first using subTri (subTri does not
% alter input geometry). This method does not alter the geometry but
% simply samples a subset of the input points. So the output point set (or
% seeds) are all points part of the original input geometry. 
% The region closest to one of the seed points can be viewed as a Voronoi
% cell. The dual of the Voronoi tesselation is a Delaunay triangulation
% which gives the output mesh connectivity. 
% The distance marching can be very slow. One tip is to do a coarse
% resampling and then to sub-triangulate the output.

%Use distance marching method
clear optionStruct
optionStruct.toleranceLevel=0; %Tolerance for convergence
optionStruct.waitBarOn=1; %Turn on/off waitbar
[Fn,Vn,seedIndex,indSeeds,d]=remeshTriSurfDistMap(F,V,numel(indListSelect)+np,indListSelect,optionStruct); %distance based marching
[~,~,ind2]=unique(seedIndex);

%%
% Visualization

cFigure;
subplot(1,3,1); hold on;
title('Distances on a triangulated surface model','fontSize',fontSize);
hp(1)=gpatch(F,V,d,'none',1); hp(1).FaceColor='Interp';
hp(2)=plotV(V(indSeeds,:),'k.','MarkerSize',markerSize);
legend(hp,{'Mesh distances','Seed point(s)'},'Location','SouthOutSide');
axisGeom;
camlight headlight;
colormap(gca,cMapDist); colorbar;

subplot(1,3,2); hold on;
title('Seed indices','fontSize',fontSize);
hp(1)=gpatch(F,V,ind2,'none',1); 
hp(2)=plotV(V(indSeeds,:),'k.','MarkerSize',markerSize);
legend(hp,{'Mesh seed indices','Seed point(s)'},'Location','SouthOutSide');
axisGeom;
camlight headlight;
colormap(gca,cMapIndices); 

subplot(1,3,3); hold on;
title('Resampled surface model','fontSize',fontSize);
plotV(V(indSeeds,:),'k.','MarkerSize',50);
hp(1)=gpatch(F,V,'kw','none',0.5);
hp(2)=gpatch(Fn,Vn,'gw','k',1,2);
legend(hp,{'Original mesh','Resampled mesh'},'Location','SouthOutSide');

axisGeom;
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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
