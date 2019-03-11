function A=polyarea_signed(V)

% function A=polyarea_signed(V)
%-------------------------------------------------------------------------
% 
%
% 
% Background (https://demonstrations.wolfram.com/SignedAreaOfAPolygon/):
% The formula for the area of a simple polygon can be elegantly derived
% using Green's theorem and extended to moments of the region. S. F.
% Bockman, "Generalizing the Formula for Areas of Polygons to Moments,"
% Amer. Math. Monthly, 96(2), 1989 pp. 131-132. 
%
%-------------------------------------------------------------------------
%%

E=[(1:size(V,1))' [(2:size(V,1))'; 1]]; %Edge array
X=V(:,1); %X coordinates
Y=V(:,2); %Y coordinates
XE=X(E);  %X coordinates of edge points
YE=Y(E); %Y coordinates of edge points
A=sum(0.5*((XE(:,1).*YE(:,2))-(XE(:,2).*YE(:,1)))); %Signed area

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
