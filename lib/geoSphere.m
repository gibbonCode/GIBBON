function [varargout]=geoSphere(varargin)

% function [F,V,Vs]=geoSphere(n,r,solidType)
% ------------------------------------------------------------------------
% This function generates an approximately geodesic triangulation of a
% sphere with radius r. The initial mesh is based on the icosahedron which
% is subdivided into triangles n times (see |subtri| function). The output
% is the triangle faces F, the Cartesian vertex coordinates V and the
% spherical vertex coordinates Vs. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 08/04/2013 Updated for GIBBON
% 12/08/2014 Added custom solidType and updated using varargin, varargout
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1 % Just n
        n=varargin{1}; %Number of refinement steps
        r=1; %Radius 
        solidType=4; %Icosahedron
    case 2
        n=varargin{1}; %Number of refinement steps
        r=varargin{2}; %Radius
        solidType=4; %Icosahedron
    case 3
        n=varargin{1}; %Number of refinement steps
        r=varargin{2}; %Radius
        solidType=varargin{3}; %Icosahedron
    otherwise
        error('False number of input arguments');
end

%%

if solidType==5
    error('Tesselation based on dodecahedron not implemented yet, use icosahedron for instance');
end

%Get initial icosahedron
[V,F]=platonic_solid(solidType,r); 

%Treat cube case 
if size(F,2)==4 
    [F,V]=quad2tri(F,V,'x');
end

%Treat dodecahedron case 
if size(F,2)==5 
    
end

% Sub-triangulate the icosahedron for geodesic sphere triangulation
if n>0 %If refinement is requested
    for q=1:n %iteratively refine triangulation and push back radii to be r        
        [F,V]=subtri(F,V,1); %Sub-triangulate      
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); %Convert to spherical coordinates
        [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,r.*ones(size(R)));  %Push back radii
    end
end

%%
varargout{1}=F;
varargout{2}=V;

if nargout==3 %Also create spherical coordinates   
    [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); 
    Vs=[T(:) P(:) R(:)];
    varargout{3}=Vs;    
end
