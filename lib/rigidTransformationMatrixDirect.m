function [varargout]=rigidTransformationMatrixDirect(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2 
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

[Q]=kabschRotationMatrix(V1,V2);

V1_m=mean(V1,1);
V2_m=mean(V2,1);

T1=eye(4,4);
T1(1:3,end)=V2_m(:);

R=eye(4,4);
R(1:3,1:3)=Q;

T2=eye(4,4);
T2(1:3,end)=-V1_m;

T=T1*R*T2;

varargout{1}=T;
varargout{2}=R(1:3,1:3);

