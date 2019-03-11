function [varargout]=groupVertices(varargin)

% function [groupIndexVertices,groupIndexFaces]=groupVertices(F,V,waitBarOption)
%-------------------------------------------------------------------------
%
%
%
% Change log: 
% 2018/08/31 Created
%
%-------------------------------------------------------------------------
%% Parse input

switch nargin
    case 1
        F=varargin{1};
        V=[];
        waitBarOption=0;
    case 2
        F=varargin{1};
        V=varargin{2};
        waitBarOption=0;
    case 3
        F=varargin{1};
        V=varargin{2};
        waitBarOption=varargin{3};
end
waitBarOption=waitBarOption>0;

if isempty(V) %If not provided assume max in F is appropriate
    numVertices=max(F(:));
else
    numVertices=size(V,1);
end

%% Get vertex-vertex connectivity

% C=patchConnectivity(F,V);
% vertexVertexConnectivity=C.vertex.vertex;
E=patchEdges(F,0); %The non-unique edge set
E_sort=sort(E,2); %Sorted in column dir so 1 2 looks the same as 2 1
indEdges=sub2indn(numVertices*ones(1,2),E_sort); %Create "virtual" indices
[~,ind1,~]=unique(indEdges); %Get indices for unique edges
E_uni=E(ind1,:); %Get unique edges
EV=[E_uni;fliplr(E_uni)];
vertexVertexConnectivity=sparse(EV(:,1),EV(:,2),EV(:,2),numVertices,numVertices);
vertexVertexConnectivity=sort(vertexVertexConnectivity,2,'descend');
[~,J,~] = find(vertexVertexConnectivity);
vertexVertexConnectivity=full(vertexVertexConnectivity(:,1:max(J)));

%% Grouping vertices

if waitBarOption==1
    hw = waitbar(0,'Grouping vertices...');
end

groupIndexVertices=zeros(numVertices,1);
groupCounter=1; %The group counter
allDone=0;
while allDone==0
    groupFinished=0;
    vertexIndicesGroup=find(groupIndexVertices==0,1); %Initialize group
    numGroupPrevious=[]; %Initialize empty group size
    while groupFinished==0
        %Grow group
        vertexIndicesGroup=vertexVertexConnectivity(vertexIndicesGroup,:);
        vertexIndicesGroup=vertexIndicesGroup(vertexIndicesGroup>0);
        vertexIndicesGroup=unique(vertexIndicesGroup);
        numGroupNow=numel(vertexIndicesGroup); %Get current group size
        if numGroupNow==numGroupPrevious
            groupFinished=1;
        end
        numGroupPrevious=numGroupNow; %Store previous group size
        if waitBarOption==1
            groupIndexVertices(vertexIndicesGroup)=groupCounter;
            waitbar(nnz(groupIndexVertices>0)/numVertices,hw);
        end
    end
    if waitBarOption~=1
        groupIndexVertices(vertexIndicesGroup)=groupCounter;
    end
    groupCounter=groupCounter+1;
    if nnz(groupIndexVertices>0)==numVertices
        allDone=1;
    end
end

if waitBarOption==1
    close(hw)
end

%% Collect output
varargout{1}=groupIndexVertices;
if nargout==2
    varargout{2}=vertexToFaceMeasure(F,groupIndexVertices);
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
