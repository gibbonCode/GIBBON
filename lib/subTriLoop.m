function [varargout]=subTriLoop(varargin)

% function [Fs,Vs,Cs,CV]=subTriLoop(F,V,n,fixBoundaryOpt,logicFixedFaces)
% ------------------------------------------------------------------------
%
%
%
%
% ------------------------------------------------------------------------

%%

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        n=1;
        fixBoundaryOpt=0;
        logicFixedFaces=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=0;
        logicFixedFaces=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=varargin{4};
        logicFixedFaces=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=varargin{4};
        logicFixedFaces=varargin{5};
    case 6
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=varargin{4};
        logicFixedFaces=varargin{5};
end

%%

if nargout>2 || ~isempty(logicFixedFaces)
    C=(1:1:size(F,1))'; %Face colors or indices
    if nargout==4
        CV=zeros(size(V,1),1); %Vertex labels/colors
    end
end

if n>0
    for qIter=1:1:n
        
        M=patchConnectivity(F,V,{'ev','ef','vv'});
        
        edgeVertexMat=M.edge.vertex;
        edgeFaceMat=M.edge.face;
        vertexVertexMat=M.vertex.vertex;
        
        if ~isempty(logicFixedFaces)
            indFixFaces=find(logicFixedFaces);            
            indFixEdges=find(any(ismember(edgeFaceMat,indFixFaces),2));            
        else
            indFixEdges=[];
        end
        
        numPoints = size(V,1);
        numEdges = size(edgeVertexMat,1);
        
        % Get indices of the three edges associated with each face
        A = sparse(edgeVertexMat(:,1),edgeVertexMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
        A = max(A,A'); %Copy symmetric
        
        %Indices for A matrix
        indA_12=F(:,1)+(F(:,2)-1)*numPoints;
        indA_23=F(:,2)+(F(:,3)-1)*numPoints;
        indA_31=F(:,3)+(F(:,1)-1)*numPoints;
        
        %Get indices for vertex array
        indV_12=full(A(indA_12));
        indV_23=full(A(indA_23));
        indV_31=full(A(indA_31));
        
        %Create faces array
        Fs=[[F(:,1)  indV_12 indV_31];...
            [F(:,2)  indV_23 indV_12];...
            [F(:,3)  indV_31 indV_23];...
            [indV_12 indV_23 indV_31]];
        
        indF1=F(edgeFaceMat(:,1),:);
        if size(edgeFaceMat,2)==2
            logicValid=edgeFaceMat(:,2)>0;
            indF2=zeros(numEdges,3);
            indF2(logicValid,:)=F(edgeFaceMat(logicValid,2),:);
        else
            indF2=zeros(numEdges,3);
        end

        
        indF=[indF1 indF2];
        
        logicPick1=indF~=edgeVertexMat(:,ones(1,6)) & indF~=edgeVertexMat(:,2*ones(1,6));
        indF(~logicPick1)=0;
        indF=sort(indF,2,'descend');
        indF=indF(:,1:2);
        
        logicBoundaryEdges=indF(:,2)==0;
        
        N=sum(vertexVertexMat>0,2);
        beta = (1./N) .* (5/8 - ( 3/8 + (1/4)*(cos(2*pi./N) )).^2);
        
        %Computed sum of neighbours
        V_sum=zeros(size(V));
        logicValid=vertexVertexMat>0;
        X=nan(size(vertexVertexMat));
        for q=1:1:size(V,2)
            X(logicValid)=V(vertexVertexMat(logicValid),q);
            V_sum(:,q)=sum(X,2,'omitnan');
        end
        
        %Compute replacement input vertices
        Vv=(1-N.*beta).*V+beta.*V_sum;
        
        if any(logicBoundaryEdges)
            eb=edgeVertexMat(logicBoundaryEdges,:);
            indBoundaryVertices=unique(eb(:));           
             
            EV=[eb;fliplr(eb)]; %Non-unique edges
            vertexVertexConnectivityBoundary=sparse(EV(:,1),EV(:,2),EV(:,2),size(V,1),size(V,1));
            vertexVertexConnectivityBoundary=sort(vertexVertexConnectivityBoundary,2,'descend');
            [~,J,~] = find(vertexVertexConnectivityBoundary);
            vertexVertexConnectivityBoundary=full(vertexVertexConnectivityBoundary(:,1:max(J)));
            
            Vv(indBoundaryVertices,:)= 6/8*V(indBoundaryVertices,:) + 1/8*(V(vertexVertexConnectivityBoundary(indBoundaryVertices,1),:)+V(vertexVertexConnectivityBoundary(indBoundaryVertices,2),:));
        end
        
        %Compute the new edge vertices        
        Vn=zeros(numEdges,size(V,2));
        if fixBoundaryOpt==1
            %Use normal mid-edge nodes for boundary edges
            Vn(logicBoundaryEdges,:)=1/2*(V(edgeVertexMat(logicBoundaryEdges,1),:) + V(edgeVertexMat(logicBoundaryEdges,2),:));
            %Replace boundary nodes with original
            indBoundaryVertices=unique(edgeVertexMat(logicBoundaryEdges,:));
            Vv(indBoundaryVertices,:)=V(indBoundaryVertices,:);
        else
            Vn(logicBoundaryEdges,:)=1/2*V(edgeVertexMat(logicBoundaryEdges,1),:) + 1/2*V(edgeVertexMat(logicBoundaryEdges,2),:);
        end
        Vn(~logicBoundaryEdges,:)=3/8*V(edgeVertexMat(~logicBoundaryEdges,1),:) + 3/8*V(edgeVertexMat(~logicBoundaryEdges,2),:) + 1/8*(V(indF(~logicBoundaryEdges,1),:)+V(indF(~logicBoundaryEdges,2),:));
        
        if ~isempty(indFixEdges)
            %Use normal mid-edge nodes for constrained edges
            Vn(indFixEdges,:)=1/2*(V(edgeVertexMat(indFixEdges,1),:) + V(edgeVertexMat(indFixEdges,2),:));
            
            %Replace constrained nodes with original
            indConstrainedVertices=unique(edgeVertexMat(indFixEdges,:));
            Vv(indConstrainedVertices,:)=V(indConstrainedVertices,:);
        end
        
        %Join point sets
        Vs = [Vv; Vn];
        
        %Color handling
        if nargout>2 || ~isempty(logicFixedFaces)
            %Override face color data
            C=repmat(C,[size(Fs,1)/size(F,1),1]);
            
            %Update vertex labels
            if nargout==4
                if qIter==1
                    CV=[zeros(size(Vv,1),1); 1*ones(size(Vn,1),1); ];
                else
                    CV=[CV; qIter.*ones(size(Vn,1),1); ];
                end
            end           
        end
        
        %Override face/vertices
        F=Fs;
        V=Vs;  
        if ~isempty(logicFixedFaces)
            logicFixedFaces=logicFixedFaces(C);
        end
    end
end

%% Collect output

varargout{1}=F;
varargout{2}=V;   
if nargout>2
    varargout{3}=C;
    if nargout==4
        varargout{4}=CV;
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
