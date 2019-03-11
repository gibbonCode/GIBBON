function [CM]=fourthOrderMat(C)

ind_C_all=1:numel(C);
[I,J,K,L]=ind2sub(size(C),ind_C_all);

%Create the 9x9 matrix indices
p=3*(I-1)+K;
q=3*(J-1)+L;

%Treat posible symbolic class
switch class(C)
    case 'sym'
        CM=sym(zeros(9,9));
    otherwise
        CM=zeros(9,9);
end

%Set values
[ind_pq]=sub2ind(size(CM),p,q);
CM(ind_pq(:))=C(:);
 
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
