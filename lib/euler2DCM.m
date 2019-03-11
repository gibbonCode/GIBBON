function [varargout]=euler2DCM(E)

% ----------------------------------------------------------------------
% function [R,Ri]=euler2DCM(E)
% 
% This function generates a rotation matrices based on the Euler angles in
% E.
%
% Kevin Mattheus Moerman, 01/06/2011
% kevinmoerman@hotmail.com
% ----------------------------------------------------------------------


switch class(E)
    case 'double'
        R=zeros(3,3,size(E,1));
    case 'sym'
        R=sym(zeros(3,3,size(E,1)));
end
Ri=R;

for q=1:1:size(E,1)
    
    Rx=[1        0              0;...
        0        cos(E(q,1))  -sin(E(q,1));...
        0        sin(E(q,1))   cos(E(q,1))];
    
    Ry=[cos(E(q,2))  0        sin(E(q,2));...
            0        1        0;...
        -sin(E(q,2)) 0        cos(E(q,2))];
    
    Rz=[cos(E(q,3))  -sin(E(q,3)) 0;...
        sin(E(q,3))  cos(E(q,3))  0;...
        0        0        1];
    
    Rxyz=Rx*Ry*Rz;
    R(:,:,q)=Rxyz;
    
    Ri(:,:,q)=inv(Rxyz);
end

varargout{1}=R; 
varargout{2}=Ri; 

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
