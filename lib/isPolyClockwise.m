function [L]=isPolyClockwise(V)

mean_V=mean(V,1);
V=V-mean_V(ones(size(V,1),1),:);

%Compute angles
T = atan2(V(:,2),V(:,1));

%Unwrap
T = unwrap(T);

%Is the mean of the derivative smaller than 0?
L=mean(diff(T))<0;