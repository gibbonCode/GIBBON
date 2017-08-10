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
% 2015/04/25 Added dodecahedron solidType
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

%Get initial icosahedron
[V,F]=platonic_solid(solidType,r); 

%Treat cube case 
if size(F,2)==4 
    [F,V]=quad2tri(F,V,'x');
end

%Treat dodecahedron case 
if size(F,2)==5 
          Vm=zeros(size(F,1),size(V,2));
        for q=1:1:size(V,2)
            X=V(:,q);
            FX=X(F);
            if size(F,1)==1 %Treat special case of single face
                FX=FX';
            end
            Vm(:,q)=mean(FX,2);
        end
        %Join point sets
        Vt=[V;Vm];
        
        indVm=(size(V,1)+1):size(Vt,1);
        
        %Create faces
        Ft=[F(:,1) F(:,2) indVm(:);...
            F(:,2) F(:,3) indVm(:);...
            F(:,3) F(:,4) indVm(:);...
            F(:,4) F(:,5) indVm(:);...
            F(:,5) F(:,1) indVm(:)];
        V=Vt;
        F=Ft;
end

if solidType~=4
    [T,P,R] = cart2sph(V(:,1),V(:,2),V(:,3)); %Convert to spherical coordinates
    [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,r.*ones(size(R)));  %Push back radii
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
 
%% <-- GIBBON footer text --> 
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
