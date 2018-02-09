function [varargout]=remeshTriSurfDistMap(varargin)
% [Fd,Vd,seedIndex,indSeed]=remeshTriSurfDistMap(F,V,numSeeds,startInds,)
% Inputs: 
% F = input mesh faces
% V = input mesh vertices
% numSeeds = number of desired vertices (must be much less than length(V))
% startInds = points to start the remeshing from (optional). These points will be included in the output mesh
% W = point weights (optional)
% Outputs: 
% Fd = Output mesh faces
% Vd = Output mesh vertices
% seedIndex = indices of selected output vertices for each index vertex
% indSeed = indices of input vertices selected for the output mesh

%% PARSE INPUT

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=1; %Use first point as only start 
        W=ones(size(V,1),1);
    case 4
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4}; 
        numStarts=numel(startInds);
        W=ones(size(V,1),1);
        if numSeeds<numStarts
            warning('numSeeds<numel(startInds) assuming numSeeds should be numSeeds=numSeeds+numel(startInds)');
            numSeeds=numSeeds+numStarts;
        end            
    case 5
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};
        numStarts=numel(startInds);
        if numSeeds<numStarts
            warning('numSeeds<numel(startInds) assuming numSeeds should be numSeeds=numSeeds+numel(startInds)');
            numSeeds=numSeeds+numStarts;
        end
        W=varargin{5};
    otherwise
        error('Wrong number of input arguments!')
end

%% COMPUTE SURFACE SEEDS AND VERTEX SEED INDICES
[indSeed,~,seedIndex]=patchMarchDistMapPointSeed(F,V,startInds,numSeeds,W);

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

varargout{1} =  Fd;
varargout{2} =  Vd;
varargout{3} =  seedIndex;
varargout{4} =  indSeed;
     
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
