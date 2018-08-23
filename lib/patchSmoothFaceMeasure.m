function [C_smooth]=patchSmoothFaceMeasure(varargin)

% function [C_smooth]=patchSmoothFaceMeasure(F,V,C,smoothPar)

%% Parse input
switch nargin 
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        smoothPar=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        smoothPar=varargin{4};
end

smoothParDefault.lambda=0.5;
smoothParDefault.n=1;
smoothParDefault.faceFaceConnectivity=[];
smoothPar=structComplete(smoothPar,smoothParDefault,1);

%% Get connectivity array
if isempty(smoothPar.faceFaceConnectivity)
    [connectivityStruct]=patchConnectivity(F,V);
    faceFaceConnectivity=connectivityStruct.face.face;
end

%%

numSmoothIterations=smoothPar.n; 
lambdaSmooth=smoothPar.lambda;

%%

nDims=size(C,2); %Number of dimensions
logicValid=faceFaceConnectivity>0;
C_smooth=C;
C_smooth_step=C; 
for qIter=1:numSmoothIterations
    %Loop for all dimensions
    for qDim=1:1:nDims
        Xp=NaN(size(C,1),size(faceFaceConnectivity,2));
        Xp(logicValid)=C_smooth(faceFaceConnectivity(logicValid),qDim);
        Xp=nanmean(Xp,2);       
        C_smooth_step(:,qDim)=Xp;
    end
    C_smooth=((1-lambdaSmooth).*C_smooth)+(lambdaSmooth.*C_smooth_step);
end
