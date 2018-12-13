function [varargout]=hemiSphereMeshHalf(varargin)

% function [F,V,C]=hemiSphereMeshHalf(nRefineSteps,sphereRadius,closeOpt)
%-------------------------------------------------------------------------
% 
%
%
% 
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        nRefineSteps=0;
        sphereRadius=1;
        closeOpt=0;        
    case 1
        nRefineSteps=varargin{1};
        sphereRadius=1;
        closeOpt=0;
    case 2
        nRefineSteps=varargin{1};
        sphereRadius=varargin{2};
        closeOpt=0;
    case 3
        nRefineSteps=varargin{1};
        sphereRadius=varargin{2};
        closeOpt=varargin{3};
end
%%

switch closeOpt
    case 0
        V=eye(3,3);
        F=1:3;
        C=1;
    case 1
        V=[eye(3,3); 0 0 0];
        F=[1 2 3; 4 2 1; 4 1 3; 4 3 2;];
        C=(1:4)';
end

for q=1:1:nRefineSteps
    
    %Refine surface
    [F,V]=subtri(F,V,1); %Sub-devide triangles
    C=repmat(C,4,1); %Replicate color data
    
    %Force nodes on sphere to conform to radius    
    switch closeOpt
        case 0
            indPush=1:1:size(V,1); %All points;
        case 1
            indPush=unique(F(C==1,:));
    end
    [azimuth,elevation,~]= cart2sph(V(indPush,1),V(indPush,2),V(indPush,3));
    [V(indPush,1),V(indPush,2),V(indPush,3)]=sph2cart(azimuth,elevation,ones(numel(indPush),1));
end

%%

%Smoothen sides
if closeOpt==1
    for q=2:1:4                                
        cPar.n=25;
        cPar.RigidConstraints=unique(F(C~=q,:));
        V=patchSmooth(F,V,[],cPar);
    end
end

%Scale by radius
V=V.*sphereRadius;

%% Collect output

varargout{1}=F;
varargout{2}=V;
if nargout>2
    varargout{3}=C;
end
