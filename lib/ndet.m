function DET=ndet(N1,N2,N3)

DET= ( N1(:,1).*N2(:,2).*N3(:,3) ) - ...
     ( N1(:,1).*N3(:,2).*N2(:,3) ) + ...
     ( N2(:,1).*N3(:,2).*N1(:,3) ) - ...
     ( N2(:,1).*N1(:,2).*N3(:,3) ) + ...
     ( N3(:,1).*N1(:,2).*N2(:,3) ) - ...
     ( N3(:,1).*N2(:,2).*N1(:,3) );
 
 
 
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
