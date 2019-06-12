function [varargout]=tri2rhombi(varargin)

% function [FQ,VQ,CQ]=tri2rhombi(F,V,C)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%% Parse input
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

%%

conStruct=patchConnectivity(F,V,{'ef','ev'});
EF=conStruct.edge.face;
E=conStruct.edge.vertex;

logicBoundary=any(EF==0,2);
indF=EF(logicBoundary,1);

VF=patchCentre(F,V);
VE=patchCentre(E,V);
VC=VF;
VC(indF,:)=VE(logicBoundary,:);

e=E(~logicBoundary,:);
ef=EF(~logicBoundary,:);

VQ=[V;VC];
FQ=[e(:,1) ef(:,1)+size(V,1) e(:,2) ef(:,2)+size(V,1)];

%Derive color data for refined set
if ~isempty(C) && nargout>2
    if size(C,1)==size(F,1) %Face color data        
        CQ=(C(ef(:,1),:)+C(ef(:,2),:))./2;
    elseif size(C,1)==size(V,1) %Vertex color data
        CC=vertexToFaceMeasure(F,C);
        CQ=[C;CC];
    else 
        error('Color data should be nxq in size whereby n is the number of faces or the number of vertices');
    end
else
    CQ=[];
end

%% Collect output
varargout{1}=FQ; %Faces
varargout{2}=VQ; %Vertices
varargout{3}=CQ; %New color data for faces or vertices

