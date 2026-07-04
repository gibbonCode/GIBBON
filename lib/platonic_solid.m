function [V,F]=platonic_solid(varargin)

% function [V,F]=platonic_solid(n,r)
% ------------------------------------------------------------------------
% PLATONIC_SOLID Creates the PATCH data, the vertices (V) and faces (F) for
% a given platonic solid (according to "n" see below) with radius (r)
%
% n=1 -> Tetrahedron
% n=2 -> Cube
% n=3 -> Octahedron
% n=4 -> Icosahedron
% n=5 -> Dodecahedron
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 12/08/2014 Updated for GIBBON
%------------------------------------------------------------------------
%% Parse input

switch nargin
    case 1
        n=varargin{1};
        r=1;
    case 2
        n=varargin{1};
        r=varargin{2};
end

if isempty(r)
    r=1; 
end

%%
switch n
    case 1 % Tetrahedron
        X=[-0.5;0.5;0;0;];
        Y=[-sqrt(3)/6;  -sqrt(3)/6; sqrt(3)/3; 0];
        Z=[-0.25.*sqrt(2/3); -0.25.*sqrt(2/3); -0.25.*sqrt(2/3);  0.75.*sqrt(2/3)];
        X=X([1 2 4 3]);
        Y=Y([1 2 4 3]);
        Z=Z([1 2 4 3]);        
        F=[3,2,1;1,2,4;2,3,4;4,3,1;];
        F=fliplr(F);
    case 2 % Cube
        X=[-1;  1; 1; -1; -1;  1; 1; -1;];
        Y=[-1; -1; 1;  1; -1; -1; 1;  1;];
        Z=[-1; -1;-1; -1;  1;  1; 1;  1;];
        F=[4,3,2,1;
            1,2,6,5;
            2,3,7,6;
            3,4,8,7;
            4,1,5,8;
            5,6,7,8;];
    case 3 % Octahedron
        X=[-1;  1; 1; -1;  0;   0;];
        Y=[-1; -1; 1;  1;  0;   0;];
        Z=[ 0;   0; 0;  0; -1;  1;];
        F=[5,2,1;5,3,2;5,4,3;5,1,4;1,2,6;2,3,6;3,4,6;4,1,6;];        
    case 4 % Icosahedron
        phi=(1+sqrt(5))/2;
        X=[0;0;0;0;-1;-1;1;1;-phi;phi;phi;-phi;];
        Y=[-1;-1;1;1;-phi;phi;phi;-phi;0;0;0;0;];
        Z=[-phi;phi;phi;-phi;0;0;0;0;-1;-1;1;1;];
        F=[9,4,1;1,5,9;1,8,5;10,8,1;4,10,1;5,2,12;12,2,3;12,3,6;12,6,9;12,9,5;10,7,11;8,10,11;2,8,11;3,2,11;7,3,11;2,5,8;10,4,7;7,6,3;6,7,4;6,4,9;];
    case 5 % Dodecahedron
        phi=(1+sqrt(5))/2;
        X=[1;(1/phi);-phi;phi;-1;0;-phi;1;-1;-1;1;(1/phi);-1;0;0;-(1/phi);phi;-(1/phi);1;0;];
        Y=[1;0;-(1/phi);(1/phi);1;-phi;(1/phi);-1;1;-1;-1;0;-1;-phi;phi;0;-(1/phi);0;1;phi;];
        Z=[1;phi;0;0;-1;-(1/phi);0;1;1;1;-1;-phi;-1;(1/phi);-(1/phi);phi;0;-phi;-1;(1/phi);];
        F=[20,9,16,2,1;2,16,10,14,8;16,9,7,3,10;7,9,20,15,5;18,13,3,7,5;3,13,6,14,10;6,13,18,12,11;6,11,17,8,14;11,12,19,4,17;1,2,8,17,4;1,4,19,15,20;12,18,5,15,19;];
    otherwise
        warning('False input for n')
end

%Altering radius
[THETA,PHI,~]=cart2sph(X,Y,Z);
[X,Y,Z]=sph2cart(THETA,PHI,r.*ones(size(X)));
V=[X(:) Y(:) Z(:)];
 
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
