function D=pathLength(V)

% function D=pathLength(V)
% ------------------------------------------------------------------------
%
% This function calculates the "current" curve patch length for each of the
% points of the curve defined by the input argument V. The curve may be
% multidimensional. The output D is a vector of size [size(V,1) 1] whereby
% each entry is defined as the sum of the point-to-point (Euclidean)
% distances (i.e. the curve patch length) leading up to that point. Hence
% it is clear that D(1) is 0 and D(end) is the total curve length. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/05/2013
%------------------------------------------------------------------------

%Compute distance metric
D=zeros(size(V,1),1); %Initialise with zeros (first values stays zero)
D(2:end)=cumsum(sqrt(sum(diff(V,1,1).^2,2)),1);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
