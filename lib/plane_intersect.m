function X=plane_intersect(V1,V2,V3,N1,N2,N3)

%%

DET=ndet(N1,N2,N3);
DET(DET==0)=NaN;

%  X= (1./DET).*[(dot(V(:,1),N(:,1)).*cross(N(:,2),N(:,3))) + (dot(V(:,2),N(:,2)).*cross(N(:,3),N(:,1))) + (dot(V(:,3),N(:,3)).*cross(N(:,1),N(:,2)))]
 
X= ((1./DET)*ones(1,3)).* (...
     ((dot(V1,N1,2)*ones(1,3)).*cross(N2,N3,2)) + ...
     ((dot(V2,N2,2)*ones(1,3)).*cross(N3,N1,2)) + ...
     ((dot(V3,N3,2)*ones(1,3)).*cross(N1,N2,2))...
     );


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
