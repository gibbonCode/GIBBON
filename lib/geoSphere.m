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
% 2013/04/08 Updated for GIBBON
% 2014/08/12 Added custom solidType and updated using varargin, varargout
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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
