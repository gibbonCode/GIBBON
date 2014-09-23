function [Vm]=rigidTransformVertices(V,T,v)

Vm=V; 
[I,J,K]=cart2im(V(:,1),V(:,2),V(:,3),v); %Convert to image coordinates
IJK=[I(:) J(:) K(:) ones(size(I(:)))]; %Prepare for mapping
IJK_mapped=(T*IJK')'; %Do mapping
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(IJK_mapped(:,1),IJK_mapped(:,2),IJK_mapped(:,3),v); %Convert mapped image coordinates back to "Cartesian" coordinates
