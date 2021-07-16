function [varargout]=subQuadCatmullClark(varargin)

% function [Fs,Vs,Cs,CV]=subQuadCatmullClark(F,V,n,fixBoundaryOpt)
% ------------------------------------------------------------------------
%
%
%
%
% ------------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        n=1;
        fixBoundaryOpt=0;
    case 3
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=0;
    case 4
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        fixBoundaryOpt=varargin{4};
end

C=(1:1:size(F,1))'; %Face colors or indices

%%

if n>0
    for qIter=1:1:n
        M=patchConnectivity(F,V,{'ev','ef','vf','ve','vv'});
        
        edgeVertexMat=M.edge.vertex;
        edgeFaceMat=M.edge.face;
        
        vertexFaceMat=M.vertex.face;
        vertexEdgeMat=M.vertex.edge;
        vertexVertexMat=M.vertex.vertex;
        
        numPoints = size(V,1);
        numEdges = size(edgeVertexMat,1);
        
        % Get indices of the three edges associated with each face
        A = sparse(edgeVertexMat(:,1),edgeVertexMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
        A = max(A,A'); %Copy symmetric
        
        %Indices for A matrix
        indA_12=F(:,1)+(F(:,2)-1)*numPoints;
        indA_23=F(:,2)+(F(:,3)-1)*numPoints;
        indA_34=F(:,3)+(F(:,4)-1)*numPoints;
        indA_41=F(:,4)+(F(:,1)-1)*numPoints;
        
        %Get indices for vertex array
        indV_12=full(A(indA_12));
        indV_23=full(A(indA_23));
        indV_34=full(A(indA_34));
        indV_41=full(A(indA_41));
        
        indV_mid=(1:1:size(F,1))'+numPoints+size(edgeVertexMat,1);
        
        %Create faces array
        Fs=[[F(:,1)  indV_12 indV_mid indV_41];...
            [F(:,2)  indV_23 indV_mid indV_12];...
            [F(:,3)  indV_34 indV_mid indV_23];...
            [F(:,4)  indV_41 indV_mid indV_34]];
        
        %Face centre points
        Vm=patchCentre(F,V);
        
        %Create vertex arrays
        
        Vne=(V(edgeVertexMat(:,1),:)+ V(edgeVertexMat(:,2),:))/2;
        
        if size(edgeFaceMat,2)==1
            logicBoundaryEdges=true(size(edgeFaceMat,1),1);
        else
            logicBoundaryEdges=edgeFaceMat(:,2)==0;
        end
        
        V_F_mean=zeros(size(V));
        logicValid=vertexFaceMat>0;
        X=nan(size(vertexFaceMat));
        for q=1:1:size(V,2)
            X(logicValid)=Vm(vertexFaceMat(logicValid),q);
            V_F_mean(:,q)=mean(X,2,'omitnan');
        end
        
        V_E_mean=zeros(size(V));
        logicValid=vertexEdgeMat>0;
        X=nan(size(vertexEdgeMat));
        for q=1:1:size(V,2)
            X(logicValid)=Vne(vertexEdgeMat(logicValid),q);
            V_E_mean(:,q)=mean(X,2,'omitnan');
        end
        
        N=sum(logicValid,2);
        Vv=(V_F_mean+2*V_E_mean+(N-3).*V)./N;
        
        
        if nnz(logicBoundaryEdges)>0
            %Use normal mid-edge nodes for boundary edges
            Vn(logicBoundaryEdges,:)=Vne(logicBoundaryEdges,:);
            indBoundaryVertices=unique(edgeVertexMat(logicBoundaryEdges,:));
            if fixBoundaryOpt==1
                %Replace boundary nodes with original
                Vv(indBoundaryVertices,:)=V(indBoundaryVertices,:);
            else
                vvm=vertexVertexMat(indBoundaryVertices,:);
                logicCheck=~ismember(vvm,indBoundaryVertices);
                vvm(logicCheck)=0;
                vvm=sort(vvm,2,'descend');
                vvm=vvm(:,[1 2]);
                
                Vv(indBoundaryVertices,:)= 6/8*V(indBoundaryVertices,:) + 1/8*(V(vvm(:,1),:)+V(vvm(:,2),:));
            end
            
        end
        
        Vn(~logicBoundaryEdges,:)=(V(edgeVertexMat(~logicBoundaryEdges,1),:)+V(edgeVertexMat(~logicBoundaryEdges,2),:)+...
            Vm(edgeFaceMat(~logicBoundaryEdges,1),:) + Vm(edgeFaceMat(~logicBoundaryEdges,2),:) )/4;
        
%         Uv=(V-Vv);
%         Un=(Uv(edgeVertexMat(:,1),:) + Uv(edgeVertexMat(:,2),:))/2;
%         Um=patchCentre(F,Uv);
%         Un(~logicBoundaryEdges,:)=(Uv(edgeVertexMat(~logicBoundaryEdges,1),:)+Uv(edgeVertexMat(~logicBoundaryEdges,2),:)+...
%             Um(edgeFaceMat(~logicBoundaryEdges,1),:) + Um(edgeFaceMat(~logicBoundaryEdges,2),:) )/4;
%         Vs = [V; Vn+Un; Vm+Um]; %Join point sets

        Vs = [Vv; Vn; Vm]; %Join point sets
        
        if qIter>1
            CV=[CV; 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1);];
        else
            CV=[zeros(size(V,1),1); 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1);];
        end
        
    %Override input for looping
    C=repmat(C,[size(Fs,1)/size(F,1),1]);
    F=Fs;
    V=Vs;
    end
    
    
else
    CV=zeros(size(V,1),1);
end

varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=CV;
