function [Dt]=patchVectorTangent(F,V,D,N)

%Derive tangent contribution through dot-product with normal vectors

if isempty(N)
    [~,~,N]=patchNormal(F,V); %Get current vertex normals    
end

Dn_mag=dot(D,N,2); %Allong normal displacement magnitudes
Dn=Dn_mag(:,ones(1,3)).*N; %Normal direction displacement vectors
Dt=D-Dn; %Tranverse or tangential only displacement vectors
