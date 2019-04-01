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

