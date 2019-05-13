function [C]=patchConnectivity(varargin)

% function [C]=patchConnectivity(F,V,conTypes)
% -----------------------------------------------------------------------
% This functions creates connectivity matrices for the input patch data
% defined by the faces F and the vertices V. The output is a structure
% containing the connectivity matrices:
%
% C.vertex.vertex
% C.vertex.face
% C.vertex.edge
% 
% C.edge.face
% C.edge.vertex
% C.edge.edge
% 
% C.face.vertex
% C.face.face
% C.face.edge
% 
% Change log: 
% 2018/08/22 
% 2019/04/22 Changed variable names to be more descriptive
% 2019/04/22 Improved function performance by creating optional output
% requests for output structure. Can be controlled through input conTypes
%
% To do: 
% Check efficiency and compare to tesIND
% -----------------------------------------------------------------------

%% Parse input
switch nargin
    case 1
        F=varargin{1};
        V=[]; 
        conTypes='all';
    case 2
        F=varargin{1};
        V=varargin{2};
        conTypes='all';
    case 3
        F=varargin{1};
        V=varargin{2};
        conTypes=varargin{3};
end

if isempty(V) %No vertices provided, assuming all points are used
    numVertices=max(F(:)); %Assume all points are used in F
else
    numVertices=size(V,1);
end
numFaces=size(F,1);
numFaceVertices=size(F,2);

%% Check requested connectivity types

conTypeSet={'vv','vf','ve','ev','ef','ee','fv','ff','fe'};
if strcmp(conTypes,'all')
    logicCompute=true(size(conTypeSet));
else
    logicCompute=contains(conTypeSet,conTypes);
end

%% Edge-vertex connectivity
if any(contains(conTypeSet(logicCompute),{'ev','ef','ve','vv','ee','ff','fe'}))
    E=patchEdges(F,0); %The non-unique edge set
    E_sort=sort(E,2); %Sorted in column dir so 1 2 looks the same as 2 1
    indEdges=sub2indn(numVertices*ones(1,2),E_sort); %Create "virtual" indices
    [~,ind1,indEdges_2]=unique(indEdges); %Get indices for unique edges
    edgeVertexConnectivity=E(ind1,:); %Get unique edges
    faceEdgeConnectivity=reshape(indEdges_2,numFaceVertices,numFaces)'; %Reshape to get results
    numEdges=size(edgeVertexConnectivity,1);
    numEdgeVertices=size(edgeVertexConnectivity,2);
end

%% Edge-face connectivity
if any(contains(conTypeSet(logicCompute),{'ef','ff'}))
    ind=(1:1:numFaces)'; %Indices for all faces
    ind=ind(:,ones(1,numFaceVertices)); %Indices copied over so it is the size of F
    ind=ind(:); %Force as column
    edgeFaceConnectivity=sparse(faceEdgeConnectivity(:),ind,ind,numEdges,numFaces,numel(ind)); %Create sparse form of connectivity matrix
    edgeFaceConnectivity=sort(edgeFaceConnectivity,2,'descend'); %Sort the sparse array
    edgeFaceConnectivity=full(edgeFaceConnectivity(:,[1 2])); %Keep relevant columns, convert to full array
end

%% Vertex-face connectivity
if any(strcmp(conTypeSet(logicCompute),'vf'))
    ind=(1:1:numFaces)';
    ind=ind(:,ones(1,numFaceVertices));
    ind=ind(:);
    vertexFaceConnectivity=sparse(F(:),ind,ind,numVertices,numFaces);
    vertexFaceConnectivity=sort(vertexFaceConnectivity,2,'descend');
    [~,J,~] = find(vertexFaceConnectivity);
    vertexFaceConnectivity=full(vertexFaceConnectivity(:,1:max(J)));
end

%% Vertex-edge connectivity
if any(contains(conTypeSet(logicCompute),{'ve','ee'}))
    ind=(1:1:numEdges)';
    ind=ind(:,ones(1,numEdgeVertices));
    ind=ind(:);
    vertexEdgeConnectivity=sparse(edgeVertexConnectivity(:),ind,ind,numVertices,numEdges);
    vertexEdgeConnectivity=sort(vertexEdgeConnectivity,2,'descend');
    [~,J,~] = find(vertexEdgeConnectivity);
    vertexEdgeConnectivity=full(vertexEdgeConnectivity(:,1:max(J)));
end

%% Vertex-vertex connectivity
if any(strcmp(conTypeSet(logicCompute),'vv'))
    EV=[edgeVertexConnectivity;fliplr(edgeVertexConnectivity)]; %Non-unique edges
    vertexVertexConnectivity=sparse(EV(:,1),EV(:,2),EV(:,2),numVertices,numVertices);
    vertexVertexConnectivity=sort(vertexVertexConnectivity,2,'descend');
    [~,J,~] = find(vertexVertexConnectivity);
    vertexVertexConnectivity=full(vertexVertexConnectivity(:,1:max(J)));
end

%% Face-face connectivity 
if any(strcmp(conTypeSet(logicCompute),'ff'))
    A=edgeFaceConnectivity(faceEdgeConnectivity(:),:);
    faceFaceConnectivity=reshape(A,numFaces,numel(A)/numFaces);
    ind=(1:1:numFaces)';
    ind=ind(:,ones(1,size(faceFaceConnectivity,2)));
    ind=ind(:);
    logicValid=faceFaceConnectivity(:)>0;
    faceFaceConnectivity=faceFaceConnectivity(logicValid);
    ind=ind(logicValid);
    faceFaceConnectivity=sparse(ind,faceFaceConnectivity(:),faceFaceConnectivity(:),numFaces,numFaces,numel(ind));
    faceFaceConnectivity(inddiag(faceFaceConnectivity))=0;
    faceFaceConnectivity=sort(faceFaceConnectivity,2,'descend');
    [~,J,~] = find(faceFaceConnectivity);
    faceFaceConnectivity=full(faceFaceConnectivity(:,1:max(J)));
end

%% Edge-edge connectivity
if any(strcmp(conTypeSet(logicCompute),'ee'))
    A=vertexEdgeConnectivity(edgeVertexConnectivity(:),:);
    edgeEdgeConnectivity=reshape(A,numEdges,numel(A)./numEdges);
    ind=(1:1:numEdges)';
    ind=ind(:,ones(1,size(edgeEdgeConnectivity,2)));
    ind=ind(:);
    logicValid=edgeEdgeConnectivity(:)>0;
    ind=ind(logicValid);
    edgeEdgeConnectivity=edgeEdgeConnectivity(logicValid);
    edgeEdgeConnectivity=sparse(ind(:),edgeEdgeConnectivity(:),edgeEdgeConnectivity(:),numEdges,numEdges,numel(ind));
    edgeEdgeConnectivity(inddiag(edgeEdgeConnectivity))=0;
    edgeEdgeConnectivity=sort(edgeEdgeConnectivity,2,'descend');
    [~,J,~] = find(edgeEdgeConnectivity);
    edgeEdgeConnectivity=full(edgeEdgeConnectivity(:,1:max(J)));
end

%% Collect output in structure

if strcmp(conTypes,'all')
    % Vertex connectivity
    C.vertex.vertex= vertexVertexConnectivity;
    C.vertex.face  = vertexFaceConnectivity;
    C.vertex.edge  = vertexEdgeConnectivity;
    % Edge connectivity
    C.edge.vertex  = edgeVertexConnectivity;
    C.edge.face    = edgeFaceConnectivity;
    C.edge.edge    = edgeEdgeConnectivity;
    % Face connectivity
    C.face.vertex  = F;
    C.face.face    = faceFaceConnectivity;
    C.face.edge    = faceEdgeConnectivity;    
else
    % Vertex connectivity
    if any(strcmp(conTypeSet(logicCompute),{'vv'}))
        C.vertex.vertex= vertexVertexConnectivity;
    end    
    if any(strcmp(conTypeSet(logicCompute),{'vf'}))
        C.vertex.face  = vertexFaceConnectivity;
    end    
    if any(strcmp(conTypeSet(logicCompute),{'ve'}))
        C.vertex.edge  = vertexEdgeConnectivity;
    end
    
    % Edge connectivity
    if any(strcmp(conTypeSet(logicCompute),{'ev'}))
        C.edge.vertex  = edgeVertexConnectivity;
    end    
    if any(strcmp(conTypeSet(logicCompute),{'ef'}))
        C.edge.face    = edgeFaceConnectivity;
    end
    if any(strcmp(conTypeSet(logicCompute),{'ee'}))
        C.edge.edge    = edgeEdgeConnectivity;
    end
    
    % Face connectivity
    if any(strcmp(conTypeSet(logicCompute),{'fv'}))
        C.face.vertex  = F;
    end
    if any(strcmp(conTypeSet(logicCompute),{'ff'}))
        C.face.face    = faceFaceConnectivity;
    end
    if any(strcmp(conTypeSet(logicCompute),{'fe'}))
        C.face.edge    = faceEdgeConnectivity;
    end
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
