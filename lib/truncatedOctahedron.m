function [FC,VC,CC]=truncatedOctahedron(r)

% function [FC,VC,CC]=truncatedOctahedron(r)
% ------------------------------------------------------------------------
% This function creates the faces (F) and vertices (V) for a
% truncated octahedron. The output C contains a labelling for the faces
% where 0 denotes hexagons and 1 denotes squares. 
% 
% 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2019/02/08 Created
% ------------------------------------------------------------------------

%% Truncating an octahedron

% Creating a platonic solid which will be truncated
[V,F]=platonic_solid(3,1.5); %Get an octahedron, the 3rd platonic solid

% Truncating a platonic solid
[FC,VC,CC]=truncatePolyhedra(F,V,1/3);

%Scale radius
[a,b,c]=cart2sph(VC(:,1),VC(:,2),VC(:,3));
[VC(:,1),VC(:,2),VC(:,3)]=sph2cart(a,b,r.*ones(size(c)));

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
