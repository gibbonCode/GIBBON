function [pp,t]=cs3dPath(V,p,w)

%Equivalent to CSCVN performance if p=1 and w=ones(size(V,1),1)

dt = sum((diff(V).^2).'); %Point spacing measure
t = cumsum([0,dt.^(1/4)]); %Curve length measure
% dt = sqrt(sum((diff(V).^2).')); %Point spacing measure
% t = cumsum([0,dt]); %Curve length measure
pp = csaps(t,V',p,[],w); %Smoothened ppform