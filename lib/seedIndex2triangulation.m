function [varargout]=seedIndex2triangulation(F,V,seedIndex)

% function [Fd,Vd,indSeed]=seedIndex2triangulation(F,V,seedIndex)
% -----------------------------------------------------------------------
%
%
%
% -----------------------------------------------------------------------

%%
if size(F,2)==4 %Quad faces
   F=[F(:,[1 2 3]);F(:,[3 4 1]);]; 
end

%%
indSeed=unique(seedIndex);
numPoints=size(V,1);
numSeeds=numel(indSeed);

%% Create seedIndex-face matrix

seedIndexFaces = seedIndex(F); %Face origin indices
seedIndexFacesSort = sort(seedIndexFaces,2); %Sort to prepare for removal of double entries
[~,ind1,~] = unique(seedIndexFacesSort,'rows'); %Find indices to remove doubles
seedIndexFaces=seedIndexFaces(ind1,:); %Remove doubles but maintain face vertex order

%% Compose Delaunay triangulation

%Find Delaunay corners where faces have 3 seed indices associated with
%their vertices i.e. where 3 voronoi cells meet
logicDelaunayCorner = seedIndexFaces(:,1)~=seedIndexFaces(:,2) & seedIndexFaces(:,2)~=seedIndexFaces(:,3) & seedIndexFaces(:,1)~=seedIndexFaces(:,3);
seedIndexFacesDelaunay=seedIndexFaces(logicDelaunayCorner,:); %Subset of Delaunay

%Check if all points are used
if ~all(ismember(indSeed,seedIndexFacesDelaunay))
    warning('Not all seed points included in Delaunay triangulation. Potential concentric Voronoi cells.')
    cleanUpNodes=true(1,1);
else
    cleanUpNodes=false(1,1);
end

%Compose Delaunay triangulation
indDelaunayFix=nan(numPoints,1);
indDelaunayFix(indSeed)=1:numSeeds;
Fd=indDelaunayFix(seedIndexFacesDelaunay); %Delaunay triangles
Vd=V(indSeed,:);  %Delaunay vertices

%Clean triangulation of potential unused points
if cleanUpNodes
    [Fd,Vd,~,~,indUni]=patchCleanUnused(Fd,Vd); %Removed unused points
    indSeed=indSeed(indUni); %Get the indices of the used seeds
end

%% Collect output
varargout{1}=Fd;
varargout{2}=Vd;
varargout{3}=indSeed;

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
