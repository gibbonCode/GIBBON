function [F,V]=quadSphere(varargin)

switch nargin
    case 1
        n=varargin{1};
        r=1; 
        oriOpt=2;
    case 2
        n=varargin{1};
        r=varargin{2};
        oriOpt=2;
    case 3
        n=varargin{1};
        r=varargin{2};
        oriOpt=varargin{3};
end
   
%% Get initial solid
% oriOpt=1 -> Tetrahedron
% oriOpt=2 -> Cube
% oriOpt=3 -> Octahedron
% oriOpt=4 -> Icosahedron
% oriOpt=5 -> Rhombic dodecahedron

switch oriOpt
    case 1 %Tetrahedron
        [V,F]=platonic_solid(1,r); 
        
        [F,V]=tri2quad(F,V);
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
    case 2 %Cube
        [V,F]=platonic_solid(2,r);
    case 3 %Octahedron
        [V,F]=platonic_solid(3,r);
        
        [F,V]=tri2quad(F,V);
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
    case 4 %Icosahedron
        [V,F]=platonic_solid(4,r);
        
        [F,V]=tri2quad(F,V);
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
    case 5 %Rhombic dodecahedron
        [F,V]=rhombicDodecahedron(r);

        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
    case 6 %geoSphere triangulation
        [F,V,~]=geoSphere(n,r);
        [F,V]=tri2quad(F,V);
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
end

%% Subquadrangulate
if oriOpt~=6
    for q=1:1:n
        [F,V]=subQuad(F,V,1);
        [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3));
        [V(:,1),V(:,2),V(:,3)]=sph2cart(T,P,ones(size(R)).*r);
    end
end
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
