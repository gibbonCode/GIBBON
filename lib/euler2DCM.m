function [varargout]=euler2DCM(E)

% function [Q,Qi]=euler2DCM(E)
% ----------------------------------------------------------------------
% 
% This function uses the input Euler angle set |a|, a 1x3 vector, to
% compute a rotation tensor |Q|, also known as a direction cosine matrix
% (DCM). 
% See also DCM2euler
%
% Change log: 
% 2011/06/01 Created
% 2019/06/21 Update variable naming and documentation
% 2019/06/21 Use transpose instead of inverse operation
% 2019/06/21 Added assumption angles are real for symbolic input
%  
% ----------------------------------------------------------------------

%% 

switch class(E)
    case 'double'
        symbolicOpt=false(1,1);
        Q=zeros(3,3,size(E,1));
    case 'sym'
        symbolicOpt=true(1,1);
        assume(E,'real'); %Assume angles are real
        assume(E>=0 & E<2*pi); %Assume range 0-2*pi
        Q=sym(zeros(3,3,size(E,1)));
end
Qi=Q;

for q=1:1:size(E,1)
    
    Qx=[1        0              0;...
        0        cos(E(q,1))  -sin(E(q,1));...
        0        sin(E(q,1))   cos(E(q,1))];
    
    Qy=[cos(E(q,2))  0        sin(E(q,2));...
            0        1        0;...
        -sin(E(q,2)) 0        cos(E(q,2))];
    
    Qz=[cos(E(q,3))  -sin(E(q,3)) 0;...
        sin(E(q,3))  cos(E(q,3))  0;...
        0        0        1];
    
    Rxyz=Qx*Qy*Qz;
    Q(:,:,q)=Rxyz;

    Qi(:,:,q)=Rxyz'; %Transpose to get inverse

end

varargout{1}=Q; 
varargout{2}=Qi; 

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
