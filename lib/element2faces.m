function [F,C]=element2faces(E,C)

%Converts elements to faces enabling PATCH based visualisation colordata is
%copied for each face. Double faces for bordering elements may occur. 
%
%
%19/01/2012, Kevin Mattheus Moerman

switch size(E,2)
    case 4 %4 node tetrahedral elements
        F=[E(:,[1 2 3]);... 
           E(:,[1 2 4]);... 
           E(:,[2 3 4]);... 
           E(:,[3 1 4])]; 
       C=repmat(C(:),4,1);
    case 8 %8 node hexahedral elements
        F=[E(:,[1 2 3 4]);... %top
           E(:,[5 6 7 8]);... %bottom
           E(:,[1 2 6 5]);... %side 1
           E(:,[3 4 8 7]);... %side 2
           E(:,[2 3 7 6]);... %front
           E(:,[1 4 8 5]);]; %back       
       C=repmat(C(:),6,1);
    otherwise
        error('MATLAB:ELEMENT2FACES: size(E,1) not consistent with tetrahedral or hexahedral element');
end

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
