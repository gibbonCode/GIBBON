function [Vm]=patchCentre(F,V)

Vm=zeros(size(F,1),3);
for q=1:1:size(V,2)
    X=V(:,q);
    Vm(:,q)=mean(X(F),2);
end