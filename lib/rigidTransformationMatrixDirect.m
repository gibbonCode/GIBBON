function [T]=rigidTransformationMatrixDirect(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2; 
    V1(:,3)=0; 
end

if size(V2,2)==2; 
    V2(:,3)=0; 
end

%% Use Kabsch algorithm

%Force input to 3D
if size(V1,2)==2; 
    V1(:,3)=0; 
end

if size(V2,2)==2; 
    V2(:,3)=0; 
end

A=V1'*V2;

[U,~,V] = svd(A);

q=sign(det(U'*V));
Q=eye(3,3);
Q(3,3)=q; 
R=V*Q*U';

T=eye(4,4);
T(1:3,1:3)=R; 
T(1:3,4)=mean(V2)-mean(V1);

