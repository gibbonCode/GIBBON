function [Fd,Vd,seedIndex]=remeshTriSurfMarchCount(varargin)


switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=1; %Use first point as only start 
    case 4
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4}; 
        numStarts=numel(startInds);
        if numSeeds<numStarts
            warning('numSeeds<numel(startInds) assuming numSeeds should be numSeeds=numSeeds+numel(startInds)');
            numSeeds=numSeeds+numStarts;
        end            
    otherwise
        error('Wrong number of input arguments!')
end


numPoints=size(V,1);
[indSeed,~,seedIndex]=patchMarchCountPointSeed(F,V,startInds,numSeeds);

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
