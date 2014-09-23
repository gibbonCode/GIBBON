function [N,Vn]=trinorm(F,V)

%N.B. if F does not describe triangles this functions uses first three
%vertices as a triangular description


%Getting triangle surface normal (cross product of two edge vectors)
vec1=[V(F(:,2),1)-V(F(:,1),1)  V(F(:,2),2)-V(F(:,1),2)  V(F(:,2),3)-V(F(:,1),3)];
vec2=[V(F(:,3),1)-V(F(:,1),1)  V(F(:,3),2)-V(F(:,1),2)  V(F(:,3),3)-V(F(:,1),3)];
N=cross(vec1,vec2,2);

%Normalizing vector length
N=N./(sqrt(sum(N.^2,2))*ones(1,size(N,2)));

%Midface coordinates for normal vectors (mean of each face)
X=V(:,1); Y=V(:,2); Z=V(:,3);
Vn=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];

end

