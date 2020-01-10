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
        Xp=gnanmean(Xp,2);       
        C_smooth_step(:,qDim)=Xp;
    end
    C_smooth=((1-lambdaSmooth).*C_smooth)+(lambdaSmooth.*C_smooth_step);
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
