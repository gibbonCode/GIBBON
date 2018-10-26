function [varargout]=edgeNormal(F,V)

% function [NE,VE]=edgeNormal(F,V)
%-------------------------------------------------------------------------
% 
% To do: Expand to add corner normals
%-------------------------------------------------------------------------

%%

%Create edge-face indices
CE=(1:size(F,1))'; %Initialize face index set
CE=[CE(:,ones(size(F,2),1))]'; %Replicate
CE=CE(:); %Force as column

%Get edge description
E=patchEdges(F,0); 

%Get edge vectors
V_edge=vecnormalize(V(E(:,2),:)-V(E(:,1),:)); %The normalized edge vectors
N_face=patchNormal(F,V); %Face normals
N_face=N_face(CE,:); %replicate using edge-face indices
NE=vecnormalize(cross(V_edge,N_face,2)); %Calculate edge normals based on cross product

%%
varargout{1}=NE;
if nargout>1
    %Compute central edge coordinates if requested
    VE=patchCentre(E,V);
    varargout{2}=VE;
end

%%