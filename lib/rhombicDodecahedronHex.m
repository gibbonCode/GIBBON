function [E,V]=rhombicDodecahedronHex(r)

%% Get rhombic dodecahedron

[F,V]=rhombicDodecahedron(r);

%% Add centre point

V(end+1,:)=zeros(1,3);

%% Define the elements

E1=[fliplr(F(1,:)) fliplr([size(V,1) F(5,3) F(5,4) F(8,1)]) ];
E2=[fliplr(F(2,:))  fliplr([F(9,2) F(10,3) size(V,1) F(9,1)]) ];
E3=[fliplr(F(3,:)) fliplr([size(V,1)  F(7,3) F(6,4) F(6,1)]) ];
E4=[fliplr((F(4,:))) fliplr([F(12,2) F(12,3) size(V,1) F(11,1)])  ];
E=[E1;E2;E3;E4];

 
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
