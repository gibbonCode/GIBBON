function [Q]=kabschRotationMatrix(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2; 
    V1(:,3)=0; 
end

if size(V2,2)==2; 
    V2(:,3)=0; 
end

%% 

%Force input to 3D
if size(V1,2)==2; 
    V1(:,3)=0; 
end

if size(V2,2)==2; 
    V2(:,3)=0; 
end

%Centering on the mean
V1_m=mean(V1,1);
V1=V1-V1_m(ones(size(V1,1),1),:);

V2_m=mean(V2,1);
V2=V2-V2_m(ones(size(V2,1),1),:);

%% Use Kabsch algorithm

A=V1'*V2;

[U,~,V] = svd(A);

q=sign(det(U'*V));
Q=eye(3,3);
Q(3,3)=q; 
Q=V*U';

