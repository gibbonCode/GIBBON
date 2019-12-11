function [U_min,U_max,C_min,C_max,C_mean,C_gauss] = patchCurvature(F,V)

% function [U_min,U_max,C_min,C_max,C_mean,C_gauss] = patchCurvature(F,V)
% ------------------------------------------------------------------------
%
% 
% 
% ------------------------------------------------------------------------

%%
% Get connectivity matrices
connecticityStruct=patchConnectivity(F,V);
edgeVertexConnectivity=connecticityStruct.edge.vertex;
edgeFaceConnectivity=connecticityStruct.edge.face;
faceEdgeConnectivity=connecticityStruct.face.edge;

numFaces = size(F,1); %Number of faces
numEdges=size(edgeVertexConnectivity,1); %Number of edges
N=patchNormal(F,V); % Face normals

% Create edge vectors
VE = V(edgeVertexConnectivity(:,1),:) - V(edgeVertexConnectivity(:,2),:); %Edge vector
edgeVectorLength = sqrt(sum(VE.^2,2)); %Length of edge vector
VE=VE./edgeVectorLength(:,ones(1,3)); %Normalized edge vector
edgeVectorLength=edgeVectorLength./mean(edgeVectorLength); %Scale edge lengths to mean

% Compute the un-signed angle
beta = real(acos(dot(N(edgeFaceConnectivity(:,1),:),N(edgeFaceConnectivity(:,2),:),2)));

% Fix sign of angle
cp = cross(N(edgeFaceConnectivity(:,1),:),N(edgeFaceConnectivity(:,2),:),2);
si = sign(dot(cp,VE,2));
beta = beta .* si;

% tensors
T = zeros(3,3,numEdges); 
for qi=1:3
    for qj=1:qi
        T(qi,qj,:) = VE(:,qi).*VE(:,qj);
        T(qj,qi,:) = T(qi,qj,:);
    end
end
betaScale=beta.*edgeVectorLength;
T = T.*repmat(permute(betaScale,[3 2 1]),[3,3,1]);

%Average edge data onto faces
TF=zeros(size(T,1),size(T,2),size(faceEdgeConnectivity,1));
for q=1:1:size(F,2)
    TF=TF+T(:,:,faceEdgeConnectivity(:,q));
end
TF=TF./size(F,2);

% Do eigen decomposition
C_min = zeros(numFaces,1); %Min eigen value
C_max = zeros(numFaces,1); %Max eigen value
U_min= zeros(numFaces,3); %Min eigen vector
U_max= zeros(numFaces,3); %Max eigen vector
for k=1:numFaces %Loop over all faces    
    %Eigen decomposition    
    [u,d] = eig(TF(:,:,k)); 
    u=real(u); 
    d = real(diag(d));
    
    %Cope with zero for normal and potential negative value for min and max
    [~,indSort1] = sort(abs(d)); % sort so zero is first followed by min, max
    [~,indSort2] = sort((d(indSort1(2:3)))); %Sort without zero
    indKeep=indSort1(2:3);
    indKeep=indKeep(indSort2);            
    C_min(k) = d(indKeep(1));
    C_max(k) = d(indKeep(2));
    U_min(k,:) = u(:,indKeep(2))';
    U_max(k,:) = u(:,indKeep(1))';
end

C_mean = (C_min+C_max)/2;
C_gauss = C_min.*C_max;

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
