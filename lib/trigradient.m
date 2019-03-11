function [ux,uy]=trigradient(TRI,V,INT)

%Equivalent to [ux,uy]=pdegrad(V',TRI',INT);


% Point indices
IND_1=TRI(:,1); IND_2=TRI(:,2); IND_3=TRI(:,3);

% Triangle sides
r23x=V(IND_3,1)-V(IND_2,1); 
r31x=V(IND_1,1)-V(IND_3,1);
r12x=V(IND_2,1)-V(IND_1,1);

r23y=V(IND_3,2)-V(IND_2,2);
r31y=V(IND_1,2)-V(IND_3,2);
r12y=V(IND_2,2)-V(IND_1,2);

ar=abs(r31x.*r23y-r31y.*r23x)/2;

g1x=-0.5*r23y./ar;
g2x=-0.5*r31y./ar;
g3x=-0.5*r12y./ar;

g1y=0.5*r23x./ar;
g2y=0.5*r31x./ar;
g3y=0.5*r12x./ar;

ux=INT(IND_1,:).*g1x + INT(IND_2,:).*g2x + INT(IND_3,:).*g3x ;
uy=INT(IND_1,:).*g1y + INT(IND_2,:).*g2y + INT(IND_3,:).*g3y;

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
