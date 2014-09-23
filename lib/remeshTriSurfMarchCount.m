function [Fd,Vd,seedIndex]=remeshTriSurfMarchCount(F,V,numSeeds)

numPoints=size(V,1);
[indSeed,~,seedIndex]=patchMarchCountPointSeed(F,V,1,numSeeds);

seedIndexFaces = seedIndex(F); %Face origin indices
seedIndexFacesSort = sort(seedIndexFaces,2); %Sort to prepare for removal of double entries
[~,ind1,~] = unique(seedIndexFacesSort,'rows'); 
seedIndexFaces=seedIndexFaces(ind1,:);%Remove doubles

%Find Voronoi corners where faces have 3 seed indices associated with their
%vertices
logicVoronoiCorner = seedIndexFaces(:,1)~=seedIndexFaces(:,2) & seedIndexFaces(:,2)~=seedIndexFaces(:,3) & seedIndexFaces(:,1)~=seedIndexFaces(:,3);
seedIndexFacesVoronoi=seedIndexFaces(logicVoronoiCorner,:);

indDelaunayFix=nan(numPoints,1);
indDelaunayFix(indSeed)=1:numSeeds;

%Compose Delaunay triangulations
Fd=indDelaunayFix(seedIndexFacesVoronoi);
Vd=V(indSeed,:);
