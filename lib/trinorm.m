function [N,Vn]=trinorm(F,V)

%N.B. if F does not describe triangles this functions uses first three
%vertices as a triangular description


%Getting triangle surface normal (cross product of two edge vectors)
vec1=[V(F(:,2),1)-V(F(:,1),1)  V(F(:,2),2)-V(F(:,1),2)  V(F(:,2),3)-V(F(:,1),3)];
vec2=[V(F(:,3),1)-V(F(:,1),1)  V(F(:,3),2)-V(F(:,1),2)  V(F(:,3),3)-V(F(:,1),3)];
N=cross(vec1,vec2,2);

%Normalizing vector length
N=N./(sqrt(sum(N.^2,2))*ones(1,size(N,2)));

%Midface coordinates for normal vectors (mean of each face)
X=V(:,1); Y=V(:,2); Z=V(:,3);
Vn=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];

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
