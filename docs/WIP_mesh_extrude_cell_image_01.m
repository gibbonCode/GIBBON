clear; close all; clc; 

%%
% Author: Kevin M. Moerman, kevin.moerman@nuigalway.ie
% Date: 2021/03/28

%% Plotting parameters

fontSize=25; 

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
loadPath=fullfile(defaultFolder,'data','PICT');
loadName='vor.jpg';

voxelSize=0.05; %The real voxel size of the image 
pointSpacing=voxelSize*18; %Point spacing to use for the mesh
layerThickness=pointSpacing*15; %Layer thickness for the mesh

t=100; %Threshold for contours
interpPar='pchip'; %Interpolation method for contour resampling
numDigitsMerge=6; %Merging of nodes similar to nth decimal place

groupingOptionStruct.outputType='label';

%% Derived parameters 
numThickness=ceil((layerThickness/pointSpacing)/3);

%% Import the structure image

m=double(importdata(fullfile(loadPath,loadName)));

%% Compute contours
% To compute contours the image is badded with either black or white
% pixels to create closed contours in MATLAB

%Creating contour set 1 (boundary cells and internal cells)
M=0*max(m(:))*ones(size(m)+2);
M(2:end-1,2:end-1)=m;
L1=M>t;
L1g=repmat(L1,1,1,3); %RGB grayscale image for plotting
[J,I]=meshgrid(1:1:size(L1,2),1:1:size(L1,1)); K=zeros(size(I));
[X,Y,Z]=im2cart(I,J,K,voxelSize);
xg=[min(X(:)) max(X(:))]; yg=[min(Y(:)) max(Y(:))];
[C1,cSiz1,cLevel1]=gcontour(X,Y,L1,0.5);
[~,indSort1]=sort(cSiz1,'descend');
C1=C1(indSort1);

%Creating contour set 2 (The connective tissue network, which coincides with internal cell boundaries)
M2=max(m(:))*ones(size(m)+2);
M2(2:end-1,2:end-1)=m;
L2=M2>t;
[C2,cSiz2,cLevel2]=gcontour(X,Y,L2,0.5);
[~,indSort2]=sort(cSiz2,'descend');
C2=C2(indSort2);

%%

cFigure; hold on; 
title('Image and raw contours');
image('XData',xg,'YData',yg,'CData',L1g);

for q=1:1:numel(C1)    
    hp1=plotV(C1{q},'g-','LineWidth',4);
end

for q=1:1:numel(C2)    
    hp2=plotV(C2{q},'r-','LineWidth',2);
end
legend([hp1 hp2],{'Cell boundary','Connective layer boundary'})
axis tight; axis equal; axis ij; 
set(gca,'FontSize',fontSize);
gdrawnow; 

%% Convert contours to vertex arrays and index vectors defining the various curve segments

[V1,indCurves1]=contour2curveSet(C1);
[V2,indCurves2]=contour2curveSet(C2);

%% Joining sets and merging the nodes 
[indCurves1,indCurves2,VC]=joinContours(indCurves1,indCurves2,V1,V2);
[indCurves1,indCurves2,VC]=mergeContours(indCurves1,indCurves2,VC,numDigitsMerge);

%% Force contours to conform to outer box
% The default contours from MATLAB do not provide a rectangular boundary as
% evident from figure one (if one zooms in on edges). This cell processses
% the outer boundary to force real corners and flat sides. 

%Corner coordinates
V_corners=[min(VC(:,1)) min(VC(:,2)); min(VC(:,1)) max(VC(:,2)); max(VC(:,1)) max(VC(:,2)); max(VC(:,1)) min(VC(:,2))];

%Use distance to find nearest corners
[d,indCorners]=minDist(V_corners,VC);

%Force points nearest to corners to be corners
VC(indCorners,:)=V_corners;

%Push edge nodes to have proper edge coordinate 
logicClose_top    = VC(:,2)>=max(VC(:,2))-voxelSize/2;
logicClose_bottom = VC(:,2)<=min(VC(:,2))+voxelSize/2;
logicClose_right  = VC(:,1)>=max(VC(:,1))-voxelSize/2;
logicClose_left   = VC(:,1)<=min(VC(:,1))+voxelSize/2;

logicAt_top    = VC(:,2)>=max(VC(:,2))-eps(0);
logicAt_bottom = VC(:,2)<=min(VC(:,2))+eps(0);
logicAt_right  = VC(:,1)>=max(VC(:,1))-eps(0);
logicAt_left   = VC(:,1)<=min(VC(:,1))+eps(0);

VC(logicClose_bottom & ~logicAt_left & ~logicAt_right,2)=min(VC(:,2));
VC(logicClose_top & ~logicAt_left & ~logicAt_right,2)=max(VC(:,2));
VC(logicClose_right & ~logicAt_top & ~logicAt_bottom,1)=max(VC(:,1));
VC(logicClose_left & ~logicAt_top & ~logicAt_bottom,1)=min(VC(:,1));

%% Resampling contours evenly based on desired point spacing
% This cell converts the curve descriptions to nx2 edge descriptions which
% are more convenient for resampling and resolving the shared nodal
% configurations. Grouping is used to split the sets into their connected
% constituents and each group is seperately resampled. Resampling behaviour
% depends on what type of curve is treated, i.e. an inner cell, a simple
% closed loop, or a boundary cell which needs its corner and edges parts
% treated separately to avoid resampling induced corner rounding. 

indBranchBoundary=indCurves2{1};
indCurves1_all=[indCurves1{:}];
indCurves2_all=[indCurves2{:}];

%Creating branch groups
E_branchBoundary=[indBranchBoundary(:) [indBranchBoundary(2:end)'; indBranchBoundary(1)]];
logicMember=ismember(E_branchBoundary,indCurves1_all);

%%
logicSideEdges=sum(logicMember,2)<2; 
E_branchBoundary_sides=E_branchBoundary(logicSideEdges,:);
E_branchBoundary_inner=E_branchBoundary(~logicSideEdges,:);

%%

[G1]=tesgroup(E_branchBoundary_sides,groupingOptionStruct);
[G2]=tesgroup(E_branchBoundary_inner,groupingOptionStruct);

E_total=[E_branchBoundary_sides;E_branchBoundary_inner];
G_total=[G1;G2+max(G1(:))];
maxIndexBranchSides=max(G1(:));
maxIndexBranches=max(G_total(:));

for q=2:1:numel(indCurves2)
    indNow=indCurves2{q};
    E_now=[indNow(:) [indNow(2:end)'; indNow(1)]];    
    E_total=[E_total;E_now];
    G_total=[G_total;(max(G_total(:))+1)*ones(size(E_now,1),1)];
end

maxIndexInner=max(G_total(:));

for q=1:1:numel(indCurves1)
    indNow=indCurves1{q};
    E_now=[indNow(:) [indNow(2:end)'; indNow(1)]]; %Current edge set    
    logicMember=ismember(E_now,indCurves2_all);  
    if ~all(logicMember(:))   
        logicKeep=sum(logicMember,2)<2; 
        E_now=E_now(logicKeep,:);
        
        logicCorner=ismember(E_now,indCorners);

        if any(logicCorner(:))
            logicNotCorner=sum(logicCorner,2)==0;
            E_not_corners=E_now(logicNotCorner,:);
            [g]=tesgroup(E_not_corners,groupingOptionStruct);
            indGroup1=unique(E_not_corners(g==min(g(:)),:));
            logicSet=any(ismember(E_now,indGroup1),2);
            E_now1=E_now(logicSet,:);
            E_now2=E_now(~logicSet,:);                        
            E_total=[E_total; E_now1; E_now2;];
            G_total=[G_total; (max(G_total(:))+1)*ones(size(E_now1,1),1); (max(G_total(:))+2)*ones(size(E_now2,1),1)];
        else
            E_total=[E_total; E_now];
            G_total=[G_total; (max(G_total(:))+1)*ones(size(E_now,1),1)];
        end
    end
end


maxIndexOuter=max(G_total(:));


cFigure; hold on; 
title('Image and contour groups');
image('XData',xg,'YData',yg,'CData',L1g);

GV=faceToVertexMeasure(E_total,VC,G_total);
gpatch(E_total,VC,GV,'interp',1,5);

axis tight; axis equal; colormap(gjet); icolorbar; axis ij; 
set(gca,'FontSize',fontSize);
gdrawnow; 

%% Process resampling
for q=1:1:max(G_total(:))
  
    logicNow=G_total==q;    
    E_now=E_total(logicNow,:); %Get current edges to resample
    
    E_total=E_total(~logicNow,:); %Remove edges from set    
    G_total=G_total(~logicNow); %Remove group from set    
    
    [E_new,V_new]=resampleEdgeSet(E_now,VC,pointSpacing,interpPar);
        
    VC=[VC; V_new]; %Add new points
    E_total=[E_total; E_new];
    G_total=[G_total; q*ones(size(E_new,1),1)];    

end

%% Remove unused points and force merging of resampled contours

[E_total,VC]=patchCleanUnused(E_total,VC);
[E_total,VC]=mergeVertices(E_total,VC);

%%

cFigure; hold on; 
title('Image and contour groups');
image('XData',xg,'YData',yg,'CData',L1g);

GV=faceToVertexMeasure(E_total,VC,G_total);
gpatch(E_total,VC,GV,'interp',1,5);

axis tight; axis equal; colormap(gjet); icolorbar; axis ij; 
set(gca,'FontSize',fontSize);
gdrawnow; 

%% Create logic arrays to capture the edges for different boundary groupings relevant for meshing

logicBoundaryCells=ismember(G_total,[maxIndexBranchSides+1:maxIndexBranches  maxIndexInner+1:maxIndexOuter]);
logicInnerCells=ismember(G_total,maxIndexBranches+1:maxIndexInner);
logicBranches=G_total<=maxIndexBranches;

%%

cFigure; hold on; 
title('Image and contour sets for meshing');
image('XData',xg,'YData',yg,'CData',L1g);
hp1=gpatch(E_total(logicBranches,:),VC,'none','b',1,6);
hp2=gpatch(E_total(logicBoundaryCells,:),VC,'none','r',1,4);
hp3=gpatch(E_total(logicInnerCells,:),VC,'none','g',1,3);
legend([hp1 hp2 hp3],{'Branch boundary','Boundary cells','Interior cells'});
axis tight; axis equal; axis ij; 
set(gca,'FontSize',fontSize);
gdrawnow; 

%% Set up mesh regions
% The single closed outer boundary is defined as a region with the cells as
% "holes" in this boundary. This defines the connective tissue mesh region.
% Next the inner cells are defined as mesh regions followed by the outer
% boundary cells. 

indNow=edgeListToCurve(E_total(logicBranches,:));
indNow=indNow(1:end-1);
Es=E_total(logicInnerCells,:); 
Gs=tesgroup(Es,groupingOptionStruct);
R=cell(1,max(Gs)+1); 
R{1,1}=VC(indNow,:);
regionSpec={};
for q=1:1:max(Gs)
    E_now=Es(Gs==q,:);
    indNow=edgeListToCurve(E_now);
    indNow=indNow(1:end-1);
    R{1,q+1}=VC(indNow,:);
    regionSpec{q+1}={VC(indNow,:);};
end
regionSpec{1}=R;
 
Es=E_total(logicBoundaryCells,:); 
Gs=tesgroup(Es,groupingOptionStruct);
for q=1:1:max(Gs)
    E_now=Es(Gs==q,:);
    indNow=edgeListToCurve(E_now);
    indNow=indNow(1:end-1);    
    regionSpec{end+1}={VC(indNow,:);};
end
    
%% Using multiRegionTriMesh2D to create the triangulated mesh

[F,V,regionLabelSurface]=multiRegionTriMesh2D(regionSpec,pointSpacing,0,0);
V(:,3)=0; %Add z coordinate
V(:,2)=-V(:,2); 
V(:,2)=V(:,2)-min(V(:,2)); 
V(:,1)=V(:,1)-min(V(:,1)); 

Eb=patchBoundary(F,V); %The outer boundary edges
Ebb=patchBoundary(F(regionLabelSurface==1,:),V); %The connective tissue boundary edges

%%

cMap=[0.5 0.5 0.5; gjet(max(regionLabelSurface)-1)]; %A custom colormap that starts grey so connective tissue stands out

cFigure; hold on; 
gpatch(F,V,regionLabelSurface);
gpatch(Eb,V,'none','k',1,4);
gpatch(Ebb,V,'none','k',1,3);
axis tight; axis equal; colormap(cMap); icolorbar; 
set(gca,'FontSize',fontSize);
drawnow; 

%% Define fibre directions for the mesh

F_con=F(regionLabelSurface==1,:); 
V_F_con=patchCentre(F_con,V);

logicCon=~all(ismember(Ebb,Eb),2);
V_Ebb=patchCentre(Ebb(logicCon,:),V);
N_Ebb=vecnormalize(V(Ebb(logicCon,1),:)-V(Ebb(logicCon,2),:));

[~,indMin]=minDist(V_F_con,V_Ebb);
N_con1=N_Ebb(indMin,:);

% nz=[0 0 1];
% N_con2=cross(N_con1,nz(ones(size(N_con1,1),1),:),2);

%%

cFigure; hold on; 
gpatch(F,V,regionLabelSurface,'none',0.5);
% gpatch(Eb,V,'none','k',1,4);
gpatch(Ebb,V,'none','k',1,3);
quiverVec(V_F_con,N_con1,pointSpacing,'k'); % plotV(V_F_con,'k.');
% quiverVec(V_F_con,N_con2,pointSpacing,'r'); % plotV(V_F_con,'k.');
axis tight; axis equal; colormap(cMap); 
set(gca,'FontSize',fontSize);
drawnow; 

%% Thicken the mesh to 3D volumetric elements
[Ep,Vp,Fp1,Fp2]=patchThick(F,V,1,layerThickness,numThickness);

%Copy region label so it holds for solid mesh
regionLabelElements=repmat(regionLabelSurface,numThickness,1); 

%Get faces for the elements for visualization
[Fp,CFp]=element2patch(Ep,regionLabelElements,'penta6');

%%

cFigure; hold on; 
gpatch(Fp,Vp,CFp,'k',0.5);

axisGeom; set(gca,'FontSize',fontSize);
colormap(cMap); icolorbar; 
drawnow; 

%% Define fibre directions for the solid mesh

V_Ep_con=patchCentre(Ep(regionLabelElements==1,:),Vp);
Np_con1=repmat(N_con1,numThickness,1);

V_Ep_cell=patchCentre(Ep(regionLabelElements~=1,:),Vp);
Np_cell=repmat([0 0 1],size(V_Ep_cell,1),1);

%%

cFigure; hold on; 
gpatch(Fp,Vp,CFp,'none',0.05);
quiverVec(V_Ep_con,Np_con1,pointSpacing,'k'); %plotV(V_Ep_con,'k.');
quiverVec(V_Ep_cell,Np_cell,pointSpacing,'r'); %plotV(V_Ep_con,'k.');
colormap(cMap); %icolorbar; 
axisGeom; set(gca,'FontSize',fontSize);
colormap(cMap); icolorbar; 
drawnow; 

%%

function [VC,indCurves]=contour2curveSet(C)

indCurves=cell(size(C));
VC=[];
indOffset=0;
for q=1:1:numel(C)
    v=C{q}; %Current contour coordinates
    n=size(v,1);
    indCurves{q}=(1:1:n)+indOffset;    
    VC=[VC; v]; 
    indOffset=indOffset+n;
end

end

function [indCurves1,indCurves2,V]=joinContours(indCurves1,indCurves2,V1,V2)

%Join node sets
V=[V1;V2]; 

%Shift indices of second set to correct for union
for q=1:1:numel(indCurves2)
    indCurves2{q}=indCurves2{q}+size(V1,1);    
end
 
end

function [indCurves1,indCurves2,V]=mergeContours(indCurves1,indCurves2,V,numDigitsMerge)

%Merge nodes
[~,indKeep,indFix]=unique(pround(V,numDigitsMerge),'rows');

%Keep unique merged nodes
V=V(indKeep,:);

%Fix indices after merging
for q=1:1:numel(indCurves1)
    indCurves1{q}=indFix(indCurves1{q})';    
end

for q=1:1:numel(indCurves2)
    indCurves2{q}=indFix(indCurves2{q})';    
end

end

function [E_new,V_new]=resampleEdgeSet(E_now,VC,pointSpacing,interpPar)

indNow=edgeListToCurve(E_now); %Curve indices for resampling
if indNow(1)==indNow(end)
    indNow=indNow(1:end-1);
    V_new=evenlySpaceCurve(VC(indNow,:),pointSpacing,interpPar,1); %Resample
    E_new=[(1:1:size(V_new,1))' [(2:1:size(V_new,1))'; 1]]+size(VC,1);
else
    V_new=evenlySpaceCurve(VC(indNow,:),pointSpacing,interpPar,0); %Resample
    E_new=[(1:1:size(V_new,1)-1)' (2:1:size(V_new,1))']+size(VC,1);
end

end


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
