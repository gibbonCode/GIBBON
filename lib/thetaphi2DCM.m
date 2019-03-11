function [R,Ri]=thetaphi2DCM(theta,phi)

Ry=[cos(phi_f) 0 sin(phi_f);...
    0 1 0;...
    -sin(phi_f) 0 cos(phi_f)];
Rz=[cos(theta_f) -sin(theta_f) 0;...
    sin(theta_f) cos(theta_f) 0;...
    0 0 1];

R=Rx*Ry*Rz;
Ri=inv(Rf);


 
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
