function D=point2TriSurfDist(F1,V1,V2)

%Get face normals
[N1,~]=trinorm(F1,V1);

%Compute face centres
X=V1(:,1); Y=V1(:,2); Z=V1(:,3);
XF=X(F1); YF=Y(F1); ZF=Z(F1);

%Position vectors for face centres
v1=[mean(XF,2) mean(YF,2) mean(ZF,2)];

%Find closest triangles
Dm=dist(V2,v1');
[~,indMin]=min(Dm,[],2); %Get indices of closest

%Difference vectors between face centres and points
x2=V2-v1(indMin,:);

%Find distance to closest triangles
D=abs(dot(N1(indMin,:),x2,2));
 
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
