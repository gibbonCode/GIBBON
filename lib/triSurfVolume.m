function [VV]=triSurfVolume(F,V)

% function [VV]=triSurfVolume(F,V)
% ------------------------------------------------------------------------
%Volume derivation based on Gauss divergence theorem
%
%%% EXAMPLE
%
% r=3; %sphere radius
% n=2; %Refinements   
% [F,V,~]=geoSphere(4,r);
% 
% VVt=(4/3)*pi*r.^3; %Theoretical volume of the sphere
% [VV]=triSurfVolume(F,V); %estimate based on triangulated surface
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 26/11/2013
%------------------------------------------------------------------------
%%

[N,~]=trinorm(F,V); %Face normals
aa=tri_area(F,V); %Areas
Zm = (V(F(:,1),3)+V(F(:,2),3)+V(F(:,3),3))./3; %Mean Z
Nz = N(:,3); %Z component of normal
vv = aa.*Zm.*Nz; %Contributions
VV = sum(abs(vv)); %Total volume

