function [I J K]=sphere_index(r,Ic,Jc,Kc)

% function [I J K]=sphere_index(r)
% ------------------------------------------------------------------------
% This function generates the indices 'I',  'J' and 'K' of voxel centers
% found within a sphere of radius 'r' (in voxels). The indices are not
% 'real' indices but are 'centered around zero'. They can be used to
% generate the indices of voxels found inside a sphere around point 'n' by
% adding In,Jn,Kn to I,J and K.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 07/08/2008
% ------------------------------------------------------------------------

if nargin == 1
    Ic=0; Jc=0; Kc=0;
end

%Set-up mesh calculate radius
[X,Y,Z] = meshgrid((-round(r+1)):1:(round(r+1)));
X=X+Jc; Y=Y+Ic; Z=Z+Kc;
radius = hypot(hypot(X,Y),Z);
[I,J,K]=ind2sub(size(X),find(radius<=r)); 
IJK_middle=round(size(X)/2);
I=I-IJK_middle(1);
J=J-IJK_middle(2);
K=K-IJK_middle(3);
 
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
