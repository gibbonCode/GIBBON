function [varargout]=meshDistMarch(varargin)

% function [d,i]=meshDistMarch(F,V,indStart,optionStruct)
% -----------------------------------------------------------------------
% The meshDistMarch function can be used to compute distances on meshes.
% The distances can be used for points sampling on the mesh or for
% remeshing. 
% The function can operate on edge descriptions or face descriptions.
% Therefore for volumetric meshes (e.g. consisting of tetrahedra or
% hexahedra) appropriate face or edge data should be computed first to
% formulate the input. 
%
% Input:
% E: the edges or faces for the mesh. E.g. an nx2 edge matrix or an nxm
% face matrix (n faces, m corners per face)
% V: the vertices for the mesh
% indStart: indices for one or more points to compute distances from
% optionStruct.toleranceLevel : The tolerance level for convergence. 0 is
% the default. 
% optionStruct.numSeeds: Defines the number of seeds to generate on the
% mesh. Default is equal to startInd. 
% optionStruct.waitBarOn=0; %Turn on/off waitbar
% optionStruct.unitEdgeOn=1; %Turn on/off the use of unit edge lengths
% 
% Output: 
% d: distances (from the start/seed points) on the mesh vertices 
% seedIndex: nearest seed (or start) point indices for each vertex, forming
% a quasi-Voronoi tesselation.
% If the second output is not requested the performance is enhanced. 
% 
% Change log:
% 2018/03/27 
% -----------------------------------------------------------------------

%% Parse input

%Default option structure
defaultOptionStruct.toleranceLevel=0;
defaultOptionStruct.waitBarOn=false(1,1);
defaultOptionStruct.unitEdgeOn=false(1,1);
defaultOptionStruct.W=[];
defaultOptionStruct.Wd=[];

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        indStart=1;
        optionStruct=defaultOptionStruct;
    case 3
        F=varargin{1};
        V=varargin{2};
        indStart=varargin{3};
        optionStruct=defaultOptionStruct;
    case 4
        F=varargin{1};
        V=varargin{2};
        indStart=varargin{3};
        optionStruct=varargin{4};
end
defaultOptionStruct.numSeeds=numel(indStart);

[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

if nargout==2
    computeSeedIndex=true(1,1);
else
    computeSeedIndex=false(1,1);
end

% Get optionional inputs
toleranceLevel=optionStruct.toleranceLevel;
numSeeds=optionStruct.numSeeds;
waitBarOn=optionStruct.waitBarOn;
unitEdgeOn=optionStruct.unitEdgeOn;

W=optionStruct.W;
if isempty(W)
    W=ones(size(V,1),1);
end

Wd=optionStruct.Wd;
if ~isempty(Wd)
    V(:,end+1)=Wd; %Add additional dimension to function as distance weight
end

numStart=numel(indStart);
if numSeeds<numStart
    warning('Number of seeds is smaller than number of start points, assuming these should be additional points')
    numSeeds=numSeeds+numStart;
    optionStruct.numSeeds=numSeeds; %override input in case function is called recursively
end
numSteps=(numSeeds-numStart)+1;

%% Variables and computations outside of iterative loop
numVertices=size(V,1); %Number of vertices
if computeSeedIndex
    indAll=(1:1:numVertices)'; %Row indices for all points
end

%Compute vertex-vertex connectivity
E=patchEdges(F,0); %The non-unique edge set
E=E(E(:,1)~=E(:,2),:); %Removed collapsed edges
E_sort=sort(E,2); %Sort in column dir so 1 2 looks the same as 2 1
indEdges=sub2indn(numVertices*ones(1,2),E_sort); %Create "virtual" indices
[~,ind1,~]=unique(indEdges); %Get indices for unique edges
E_uni=E(ind1,:); %Get unique edges
EV=[E_uni;fliplr(E_uni)]; 
C=sparse(EV(:,1),EV(:,2),EV(:,2),numVertices,numVertices);
C=sort(C,2,'descend');
[~,J,~] = find(C);
C=full(C(:,1:max(J)));
siz=size(C); %Size of C
L=C>0; %Logic for valid indices in C

%Calculate edge lengths
if unitEdgeOn
    DE=nan(siz);
    DE(L)=1; %Unit edge lenghts
else
    DE=zeros(siz);
    for q=1:1:size(V,2)
        X=V(:,q);
        V_X=X(:,ones(siz(2),1));
        V_X(~L)=nan;
        VV_X=nan(siz);
        VV_X(L)=X(C(L));
        DE=DE+(VV_X-V_X).^2;
    end
    DE=sqrt(DE);
end

%Initialize distance vector
d=nan(numVertices,1);
d(indStart)=0;
DD=nan(siz);
d_previous=d;
IND_d=C(L);

%Initialize seed index vector
if computeSeedIndex
    i=nan(numVertices,1);
    i(indStart)=indStart;
    II=nan(siz);
end

%% Propagate and iteratively update distances
indSeed=nan(1,numSeeds);
indSeed(1:numStart)=indStart;

if numSteps>1      
    numLoopSteps=numSteps+1;
else
    numLoopSteps=1;
end

if waitBarOn
    hw=waitbar(1/numLoopSteps,['meshDistMarch...',num2str(round(100.*1/numLoopSteps)),'%']);
end

for q=1:1:numLoopSteps
    indStart=indSeed(~isnan(indSeed));
    while 1
        %Update distance data
        DD(L)=d(IND_d); %The distance data currently at neighbours
        D_check=DD+DE; %Distance data plus edge lenghts
        
        if computeSeedIndex
            [d,indMin]=min(D_check,[],2,'omitnan'); %Assign minimum distance
            
            %Update index data
            II(L)=i(IND_d); %The seed indices currently at neighbours
            ind=sub2ind(siz,indAll,indMin); %Linear indices for minimum used
            i=II(ind); %Get seed index at point chosen
            i(indStart)=indStart; %Override start seed indices
        else
            d=min(D_check,[],2,'omitnan'); %Assign minimum distance
        end
        d(indStart)=0; %Override starts to be zero

        %Check convergence
        if nnz(isnan(d))==nnz(isnan(d_previous))  %Check if number of NaN's changed
            if sum((d_previous-d).^2,'omitnan')<=toleranceLevel %Check if converged to within tolerance
                break %Stop while loop
            end
        end
        d_previous=d; %Store previous
    end
    if q>1 && q<=numSteps
        [~,indSeed(numStart+(q-1))]=max(d.*W,[],'omitnan');
    end
    
    if waitBarOn
        waitbar(q/numLoopSteps,hw,['meshDistMarch...',num2str(round(100.*q/numLoopSteps)),'%']);
    end    
end

if waitBarOn
    close(hw);
end

%% Store output
varargout{1}=d;
if nargout==2
    varargout{2}=i;
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
