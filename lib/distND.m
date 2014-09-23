function [D]=distND(V1,V2)

D=zeros(size(V1,1),size(V2,1));

for q=1:1:size(V1,2) %For all dimensions
    A=V1(:,q);
    B=V2(:,q);
    
    AB=A(:,ones(1,size(B,1)));
    BA=B(:,ones(1,size(A,1)))';
    
    D=D+(AB-BA).^2;
end
D=sqrt(D);