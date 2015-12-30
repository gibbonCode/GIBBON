function [T]=affineTransformationMatrixDirect(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2; 
    V1(:,3)=0; 
end

if size(V2,2)==2; 
    V2(:,3)=0; 
end

%Expand to nx4
V1_M=V1;
V1_M(:,4)=1; 
V2_M=V2;
V2_M(:,4)=1; 

%Get transformation using left devide
T=(V1_M\V2_M)';

