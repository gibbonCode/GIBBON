function [varargout]=patchBoundary(F,V)

%Get non-unique edges
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

%Get boundary indices
[indBoundary]=tesBoundary(E,V);

%Boundary edges
Eb=E(indBoundary,:);

%Output
varargout{1}=Eb; %Boundary edges
varargout{2}=E; %All edges
varargout{3}=indBoundary; %Indices for boundary edges