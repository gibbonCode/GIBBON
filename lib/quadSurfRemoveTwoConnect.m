function [varargout]=quadSurfRemoveTwoConnect(varargin)

% function [F,V,C,indFix,L,logicTwoConnect]=quadSurfRemoveTwoConnect(F,V,C)
% -----------------------------------------------------------------------
%
%
%
%
% -----------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
end

if isempty(C)
    C=ones(size(F,1),1);
end

%%
conStruct=patchConnectivity(F,V,'vf');
conVertexFace=conStruct.vertex.face;

Eb=patchBoundary(F,V); %Boundary edges
logicBoundary=false(size(V,1),1);
logicBoundary(Eb(:))=1;
logicTwoConnect=(sum(conVertexFace>0,2)<=2) & ~logicBoundary;

indVertices=find(logicTwoConnect);
Fqn=zeros(nnz(logicTwoConnect),size(F,2));
Cqn=zeros(nnz(logicTwoConnect),size(C,2));
for q=1:1:nnz(logicTwoConnect)
    indVertexNow=indVertices(q);
    indFacesNow=conVertexFace(indVertexNow,:);
    indFacesNow=indFacesNow(indFacesNow>0);
    E=patchEdges(F(indFacesNow,:));
    logicRemove=any(E==indVertexNow,2);
    E_cand=E(~logicRemove,:);
    E1=E_cand(1,:);
    Fn=[E1 E_cand(~any(ismember(E_cand,E1),2),:)];
    
    
    Fqn(q,:)=Fn;
    Cqn(q,:)=mean(C(indFacesNow,:),1);
end
indFacesRemove=conVertexFace(logicTwoConnect,:);
indFacesRemove=indFacesRemove(indFacesRemove>0);

logicFacesKeep=true(size(F,1),1);
logicFacesKeep(indFacesRemove)=0;

F=[F(logicFacesKeep,:);Fqn];
C=[C(logicFacesKeep,:);Cqn];
L=[false(nnz(logicFacesKeep),1); true(nnz(logicTwoConnect),1)];

[F,V,indFix]=patchCleanUnused(F,V);

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=indFix;
varargout{5}=L;
varargout{6}=logicTwoConnect;

