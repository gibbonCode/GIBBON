function [Et,Vt,Ct]=octa2tet(E,V)


%%

%Get octahedron element faces
C=(1:1:size(E,1))';
[F,CF]=element2patch(E,C,'octa6');

%Get centre coordinates
Vn=patchCentre(E,V);

%Append new centre coordinates
Vt=[V;Vn]; 

%Create new tetrahedral element array 
indAdd=((size(V,1)+1:size(Vt,1))'*ones(1,8));
Et=[indAdd(:) F];

%Gather color information
Ct=CF;

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
