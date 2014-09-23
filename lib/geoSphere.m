function [F,V,Vs]=geoSphere(n,r)

% function [F,V]=geoSphere(n,r)
% ------------------------------------------------------------------------
% This function generates an approximately geodesic triangulation of a
% sphere with radius r. The initial mesh is based on the icosahedron which
% is subdivided into triangles n times (see |subtri| function). The output
% is the triangle faces F, the Cartesian vertex coordinates V and the
% spherical vertex coordinates Vs. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/04/2013
%------------------------------------------------------------------------

%Get initial icosahedron
[V,F]=platonic_solid(4,r); 

% Sub-triangulate the icosahedron for geodesic sphere triangulation
if n>0 %If refinement is requested
    for q=1:n %iteratively refine triangulation and push back radii to be r        
        [F,V]=subtri(F,V,1); %Sub-triangulate      
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); %Convert to spherical coordinates
        [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,r.*ones(size(R)));  %Push back radii
    end
end
[T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); %Spherical coordinates
Vs=[T(:) P(:) R(:)];
