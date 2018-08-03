function A=polyarea_signed(V)

% function A=polyarea_signed(V)
%-------------------------------------------------------------------------
% 
%
% 
% Background (https://demonstrations.wolfram.com/SignedAreaOfAPolygon/):
% The formula for the area of a simple polygon can be elegantly derived
% using Green's theorem and extended to moments of the region. S. F.
% Bockman, "Generalizing the Formula for Areas of Polygons to Moments,"
% Amer. Math. Monthly, 96(2), 1989 pp. 131-132. 
%
%-------------------------------------------------------------------------
%%

E=[(1:size(V,1))' [(2:size(V,1))'; 1]]; %Edge array
X=V(:,1); %X coordinates
Y=V(:,2); %Y coordinates
XE=X(E);  %X coordinates of edge points
YE=Y(E); %Y coordinates of edge points
A=sum(0.5*((XE(:,1).*YE(:,2))-(XE(:,2).*YE(:,1)))); %Signed area

