function [Fd,Vd,seedIndex]=remeshTriSurfDistMapIterative(varargin)

%% PARSE INPUT

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=1; %Use first point as only start
        W=ones(size(V,1),1);
        optStruct.waitBarOn=1;
        optStruct.numIterationsMax=500;
        optStruct.toleranceLevel=1e-4;
    case 4
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};
        W=ones(size(V,1),1);
        optStruct.waitBarOn=1;
        optStruct.numIterationsMax=500;
        optStruct.toleranceLevel=1e-4;
    case 5
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};
        W=varargin{5};
        optStruct.waitBarOn=1;
        optStruct.numIterationsMax=500;
        optStruct.toleranceLevel=1e-4;
    case 6
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};
        W=varargin{5};
        optStruct=varargin{6};
    otherwise
        error('Wrong number of input arguments!')
end

numStarts=numel(startInds);
if numSeeds<numStarts
            warning('numSeeds<numel(startInds) assuming numSeeds should be numSeeds=numSeeds+numel(startInds)');
            numSeeds=numSeeds+numStarts;
end

%% COMPUTE SURFACE SEEDS AND VERTEX SEED INDICES

[indSeed,~,seedIndex]=patchMarchDistMapPointSeedIterative(F,V,startInds,numSeeds,W,optStruct);

numPoints=size(V,1);

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
