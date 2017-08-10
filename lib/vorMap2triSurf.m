function [Fd,Vd]=vorMap2triSurf(F,V,indSeed,seedIndex)

numPoints=size(V,1);
numSeeds=numel(indSeed);

seedIndexFaces = seedIndex(F); %Face origin indices
seedIndexFacesSort = sort(seedIndexFaces,2); %Sort to prepare for removal of double entries
[~,ind1,~] = unique(seedIndexFacesSort,'rows'); %Find indices to remove doubles
seedIndexFaces=seedIndexFaces(ind1,:); %Remove doubles but maintain face vertex order

%Find Voronoi corners where faces have 3 seed indices associated with their vertices
logicVoronoiCorner = seedIndexFaces(:,1)~=seedIndexFaces(:,2) & seedIndexFaces(:,2)~=seedIndexFaces(:,3) & seedIndexFaces(:,1)~=seedIndexFaces(:,3);
seedIndexFacesVoronoi=seedIndexFaces(logicVoronoiCorner,:);

%Compose Delaunay triangulation
indDelaunayFix=nan(numPoints,1);
indDelaunayFix(indSeed)=1:numSeeds;

%Prepare output
Fd=indDelaunayFix(seedIndexFacesVoronoi); %Delaunay triangular faces
Vd=V(indSeed,:);  %Delaunay vertices
 
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
