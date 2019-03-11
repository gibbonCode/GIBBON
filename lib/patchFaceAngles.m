function [betaAngle] = patchFaceAngles(F,V)

% function [betaAngle] = patchFaceAngles(F,V)
% ------------------------------------------------------------------------
%
%
% 2018/11/14 Created
% 
% ------------------------------------------------------------------------

%%
% Get connectivity matrices
connecticityStruct=patchConnectivity(F,V);
edgeVertexConnectivity=connecticityStruct.edge.vertex;
edgeFaceConnectivity=connecticityStruct.edge.face;
faceEdgeConnectivity=connecticityStruct.face.edge;
% vertexEdgeConnectivity=connecticityStruct.vertex.edge;

logicValid=all(edgeFaceConnectivity>0,2);

% Face normals
N=patchNormal(F,V); 

% Create edge vectors
VE = V(edgeVertexConnectivity(:,1),:) - V(edgeVertexConnectivity(:,2),:); %Edge vector
edgeVectorLength = sqrt(sum(VE.^2,2)); %Length of edge vector
VE=VE./edgeVectorLength(:,ones(1,3)); %Normalized edge vector


% Compute the un-signed angle
betaAngle=nan(size(edgeFaceConnectivity,1),1); 
betaAngle(logicValid) = real(acos(dot(N(edgeFaceConnectivity(logicValid,1),:),N(edgeFaceConnectivity(logicValid,2),:),2)));

% Fix sign of angle
cp = nan(size(edgeFaceConnectivity,1),3); 
cp(logicValid,:) = cross(N(edgeFaceConnectivity(logicValid,1),:),N(edgeFaceConnectivity(logicValid,2),:),2);

si = sign(dot(cp,VE,2));
% si(abs(si)<eps)=1;

betaAngle=pi+(betaAngle.*si);

betaAngle=betaAngle(faceEdgeConnectivity);

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
