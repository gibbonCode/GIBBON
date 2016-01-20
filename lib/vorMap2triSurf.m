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