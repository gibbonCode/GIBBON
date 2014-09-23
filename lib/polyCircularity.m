function C=polyCircularity(V)

%Compute area
A=patch_area(1:size(V,1),V);

%Compute circularity
Vc=[V; V(1,:)]; %Close loop
D=pathLength(Vc);
D=max(D);
C=4*pi.*(A./(D.^2));