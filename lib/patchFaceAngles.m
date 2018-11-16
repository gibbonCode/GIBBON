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