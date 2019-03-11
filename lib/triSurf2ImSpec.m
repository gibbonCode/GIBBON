function [M,G]=triSurf2ImSpec(F,V,voxelSize)

% function [M,G,bwLabels]=triSurf2Im(F,V,voxelSize)
% -----------------------------------------------------------------------
% This function converts the input triangulated surface, specified by the
% faces F and the vertices V into an image based on the voxel size
% specified in the input voxelSize. If the latter is empty then the
% voxelSize is based on the maximum edge length. If required the surface is
% resampled to subvoxel resolution using the |subtri| function. Then
% surface vertices are simply mapped to an image coordinate system and
% since they are densely sampled with respect to the voxel size form an
% enclosing boundary of voxels. The exterior, boundary and interior are
% then formulated using the bwlabeln function. The exterior, boundary and
% intertior voxels are labelled using 0's, 1's and 2's in the output image
% M. The second output G is a structure containing the fields G.voxelSize
% and G.origin which form the geometric information for the image. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/08/2013
%------------------------------------------------------------------------

%%
%Checking edge lenghts of surface
TR = triangulation(F,V);
E=edges(TR);
Ed=0;
for q=1:1:size(V,2)
    Ed=Ed+(V(E(:,1),q)-V(E(:,2),q)).^2;
end
edgeLengths=sqrt(Ed);
maxEdgeLength=max(edgeLengths(:));

if isempty(voxelSize)
    voxelSize=maxEdgeLength;    
end

%Resample surface if voxelsize is small with respect to edgelenghts
n=maxEdgeLength/voxelSize;
if n>(1-eps(1))
    n=floor(n);
    [~,V]=subtri(F,V,n);
end

%Determine surface set coordinate minima
minV=min(V,[],1);

%Determine shift so all coordinates are positive
imOrigin=-(minV-voxelSize);

%Shift points using origin
V=V+imOrigin(ones(size(V,1),1),:);

%Convert to image coordinates
V_IJK=V;
[V_IJK(:,1),V_IJK(:,2),V_IJK(:,3)]=cart2im(V(:,1),V(:,2),V(:,3),voxelSize*ones(1,3));

%Rounding image coordinates to snap to voxel
V_IJK=round(V_IJK);

%Determin image size
siz=max(V_IJK,[],1)+1;

%Get linear indices of points
indV=sub2ind(siz,V_IJK(:,1),V_IJK(:,2),V_IJK(:,3));

%Create surface boundary image
L=false(siz);
L(indV)=1;

%Create boundary, interior and exterior image  
LL = bwlabeln(~L,6); %Get labels for ~L image which will segment interior, exterior and boundary
uniqueLabels=unique(LL(:));

labelsBoundary=LL(L); %The label numbers for the boundary

indExteriorVoxel=1; %First is outside since image is at least a voxel too big on all sides
labelExterior=LL(indExteriorVoxel); %The label number for the exterior

labelsInterior=uniqueLabels(~ismember(uniqueLabels,[labelsBoundary(:); labelExterior])); %Labels for the interior (possibly multiple)

M=zeros(size(L)); %The exterior is set to 0
M(L)=1; %The boundary is set to 1
M(ismember(LL,labelsInterior))=2; %Interior is set to 2

%Storing image geometry metrics
G.voxelSize=voxelSize;
G.origin=imOrigin;
 
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
