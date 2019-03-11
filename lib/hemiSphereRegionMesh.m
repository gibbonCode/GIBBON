function [F,V,regionInd]=hemiSphereRegionMesh(inputStruct)

% function [F,V,regionInd]=hemiSphereRegionMesh(inputStruct)
% ------------------------------------------------------------------------
% This function generates an approximately geodesic triangulated mesh on a
% hemisphere whereby the mesh is grouped into triangular regions.
%
% The following inputs are required:
% hemiSphereStruct.sphereRadius -> The sphere radius
% hemiSphereStruct.nRefineRegions -> Number of refinement steps for regions
% hemiSphereStruct.nRefineMesh -> Number of refinement steps for mesh
%
% The mesh starts of as a cropped icosahedron with 12 triangular faces.
% Each refinement step introduces 4 new triangular faces per original
% triangular face.
%
% The outputs are:
% F -> the triangular faces
% V -> the vertex coordinates
% regionIndSub -> the region index numbers
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2013/04/08
% 2016/09/12 Updated documentation and integration with GIBBON
% 2018/08/19 Removed smoothing requirement
% 2018/08/19 Use hemiSphereMesh instead
% 2018/08/19 Improved input and output handling
% 2018/08/19 Added close option
%------------------------------------------------------------------------

%% Parse input

defaultInput.sphereRadius=1;
defaultInput.nRefineRegions=0;
defaultInput.nRefineMesh=1;
defaultInput.closeOption=0; 
[inputStruct]=structComplete(inputStruct,defaultInput,1);

%Sub-triangulation refinement introduces 4 triangles for each triangle per iteration  
sphereRadius=inputStruct.sphereRadius; %Sphere radius
nRefineRegions=inputStruct.nRefineRegions; %Number of refinement steps for regions
nRefineMesh=inputStruct.nRefineMesh; %Number of refinement steps for mesh
closeOption=inputStruct.closeOption; %Option to close the bottom of the hemisphere

%% Get hemi-sphere mesh

[F,V]=hemiSphereMesh(nRefineRegions,sphereRadius,closeOption);
regionInd=(1:1:size(F,1))'; %Initialize region indices as face indices

%% Sub-triangulate

for q=1:1:nRefineMesh   
    [F,V]=subtri(F,V,1); %Sub-triangulate surface    
    [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
    [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,sphereRadius*ones(size(R))); % Push back radii
    regionInd=repmat(regionInd,4,1); %Replicate color data  
end

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
